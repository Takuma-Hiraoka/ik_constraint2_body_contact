#include <ik_constraint2_body_contact/BodyContactConstraint.h>
#include <ik_constraint2/Jacobian.h>
#include <cnoid/TimeMeasure>

namespace ik_constraint2_body_contact{
  size_t getJointDOF(const cnoid::LinkPtr& joint) {
    if(joint->isRevoluteJoint() || joint->isPrismaticJoint()) return 1;
    else if(joint->isFreeJoint()) return 6;
    else return 0;
  }

  unsigned int BodyContactConstraint::convertContactPointsIdx(int x, int y, int z, int dim) {
    return (x+dim/2) + dim*(y+dim/2) + dim*dim*(z+dim/2);
  }

  void BodyContactConstraint::A_setContactPoints(std::vector<cnoid::Isometry3> contactPoints, double contactPointLength, int contactPointAreaDim) {
    this->A_contactPointLength_ = contactPointLength;
    this->A_contactPointAreaDim_ = contactPointAreaDim;
    this->A_contactPoints_.resize(this->A_contactPointAreaDim_*this->A_contactPointAreaDim_*this->A_contactPointAreaDim_);
    for (int i=0; i<contactPoints.size(); i++) {
      int x = std::floor(contactPoints[i].translation()[0]/this->A_contactPointLength_);
      int y = std::floor(contactPoints[i].translation()[1]/this->A_contactPointLength_);
      int z = std::floor(contactPoints[i].translation()[2]/this->A_contactPointLength_);
      if (x>= this->A_contactPointAreaDim_/2 || x < -this->A_contactPointAreaDim_/2 ||
          y>= this->A_contactPointAreaDim_/2 || y < -this->A_contactPointAreaDim_/2 ||
          z>= this->A_contactPointAreaDim_/2 || z < -this->A_contactPointAreaDim_/2) {
        std::cerr << "[BodyContactConstraint] contactPoint is out of contactPointArea. Check contactPointLength and contactPointAreaDim !! " << contactPoints[i].translation().transpose() << std::endl;
      }
      this->A_contactPoints_[convertContactPointsIdx(x,y,z, this->A_contactPointAreaDim_)].push_back(contactPoints[i]);
    }
  }

  void BodyContactConstraint::B_setContactPoints(std::vector<cnoid::Isometry3> contactPoints, double contactPointLength, int contactPointAreaDim) {
    this->B_contactPointLength_ = contactPointLength;
    this->B_contactPointAreaDim_ = contactPointAreaDim;
    this->B_contactPoints_.resize(this->B_contactPointAreaDim_*this->B_contactPointAreaDim_*this->B_contactPointAreaDim_);
    for (int i=0; i<contactPoints.size(); i++) {
      int x = std::floor(contactPoints[i].translation()[0]/this->B_contactPointLength_);
      int y = std::floor(contactPoints[i].translation()[1]/this->B_contactPointLength_);
      int z = std::floor(contactPoints[i].translation()[2]/this->B_contactPointLength_);
      if (x>= this->B_contactPointAreaDim_/2 || x < -this->B_contactPointAreaDim_/2 ||
          y>= this->B_contactPointAreaDim_/2 || y < -this->B_contactPointAreaDim_/2 ||
          z>= this->B_contactPointAreaDim_/2 || z < -this->B_contactPointAreaDim_/2) {
        std::cerr << "[BodyContactConstraint] contactPoint is out of contactPointArea. Check contactPointLength and contactPointAreaDim !! " << contactPoints[i].translation().transpose() << std::endl;
      }
      this->B_contactPoints_[convertContactPointsIdx(x,y,z, this->B_contactPointAreaDim_)].push_back(contactPoints[i]);
    }
  }

  void BodyContactConstraint::calcNominals(const std::vector<std::vector<cnoid::Isometry3> >& cps, // input
                                           const double& cpLength,
                                           const int& cpAreaDim,
                                           std::vector<std::vector<cnoid::Vector3> >& normals, // output
                                           std::vector<std::vector<cnoid::Vector3> >& normalJXs,
                                           std::vector<std::vector<cnoid::Vector3> >& normalJYs,
                                           std::vector<std::vector<cnoid::Vector3> >& normalJZs) {
    normals.resize(cpAreaDim*cpAreaDim*cpAreaDim);
    normalJXs.resize(cpAreaDim*cpAreaDim*cpAreaDim);
    normalJYs.resize(cpAreaDim*cpAreaDim*cpAreaDim);
    normalJZs.resize(cpAreaDim*cpAreaDim*cpAreaDim);
    double normalAngle = M_PI * 2 / 3; // 薄い板の場合の裏側は除く
    for (unsigned int i=0; i<cps.size(); i++) {
      if (cps[i].size() == 0) continue;
      // 計算する近傍のidxを求める
      int x = i % cpAreaDim - (cpAreaDim)/2;
      int y = (i % (cpAreaDim*cpAreaDim)) / cpAreaDim - (cpAreaDim)/2;
      int z = i / (cpAreaDim*cpAreaDim) - (cpAreaDim)/2;
      std::vector<unsigned int> idxs;
      int length = this->normalGradientDistance_ / cpLength;
      idxs.push_back(i);
      for (int kx=-length; kx<=length;kx++){
        for (int ky=-length; ky<=length;ky++){
          for (int kz=-length; kz<=length;kz++){
            unsigned int k = convertContactPointsIdx(x+kx,y+ky,z+kz, cpAreaDim);
            if (k<cps.size()) idxs.push_back(k);
          }
        }
      }

      for (int j=0; j<cps[i].size(); j++){
        cnoid::Vector3 normalSum = cnoid::Vector3::Zero();
        cnoid::Vector3 normalJXSum = cnoid::Vector3::Zero();
        cnoid::Vector3 normalJYSum = cnoid::Vector3::Zero();
        cnoid::Vector3 normalJZSum = cnoid::Vector3::Zero();
        int weightSum = 0;
        for (unsigned int idx=0; idx<idxs.size(); idx++){
          for (int k=0;k<cps[idxs[idx]].size(); k++) {
            double dist = (cps[i][j].translation() - cps[idxs[idx]][k].translation()).norm();
            if (dist < this->normalGradientDistance_ &&
                (normalAngle >= std::acos(std::min(1.0,(std::max(-1.0,(cps[i][j].linear()*cnoid::Vector3::UnitZ()).dot(cps[idxs[idx]][k].linear() * cnoid::Vector3::UnitZ()))))))) {
              normalSum += cps[idxs[idx]][k].linear()*cnoid::Vector3::UnitZ();
              if ((cps[idxs[idx]][k].translation() - cps[i][j].translation())[0] != 0) normalJXSum += (cps[i][j].linear()*cnoid::Vector3::UnitZ()).cross(cps[idxs[idx]][k].linear()*cnoid::Vector3::UnitZ()) / (cps[idxs[idx]][k].translation() - cps[i][j].translation())[0];
              else normalJXSum += cnoid::Vector3::Zero();
              if ((cps[idxs[idx]][k].translation() - cps[i][j].translation())[1] != 0) normalJYSum += (cps[i][j].linear()*cnoid::Vector3::UnitZ()).cross(cps[idxs[idx]][k].linear()*cnoid::Vector3::UnitZ()) / (cps[idxs[idx]][k].translation() - cps[i][j].translation())[1];
              else normalJYSum += cnoid::Vector3::Zero();
              if ((cps[idxs[idx]][k].translation() - cps[i][j].translation())[2] != 0) normalJZSum += (cps[i][j].linear()*cnoid::Vector3::UnitZ()).cross(cps[idxs[idx]][k].linear()*cnoid::Vector3::UnitZ()) / (cps[idxs[idx]][k].translation() - cps[i][j].translation())[2];
              else normalJZSum += cnoid::Vector3::Zero();
              weightSum++;
            }
          }
        }
        normals[i].push_back((normalSum / weightSum).normalized());
        normalJXs[i].push_back((normalJXSum / weightSum));
        normalJYs[i].push_back((normalJYSum / weightSum));
        normalJZs[i].push_back((normalJZSum / weightSum));
      }
    }
  }

  void BodyContactConstraint::selectNearestPoint(const std::vector<std::vector<cnoid::Isometry3> >& cps, // input
                                                 const double& cpLength,
                                                 const double& cpAreaDim,
                                                 const std::vector<std::vector<cnoid::Vector3> >& normals,
                                                 const std::vector<std::vector<cnoid::Vector3> >& normalJXs,
                                                 const std::vector<std::vector<cnoid::Vector3> >& normalJYs,
                                                 const std::vector<std::vector<cnoid::Vector3> >& normalJZs,
                                                 const cnoid::LinkPtr cp_link, // output
                                                 cnoid::Vector3& normal,
                                                 std::vector<cnoid::Vector3>& normalJacobian,
                                                 cnoid::Isometry3& localPos) {
    int x = cp_link->p()[0]/cpLength;
    int y = cp_link->p()[1]/cpLength;
    int z = cp_link->p()[2]/cpLength;
    int dx = cp_link->p()[0] - x*cpLength > cpLength/2 ? 1 : -1;
    int dy = cp_link->p()[1] - y*cpLength > cpLength/2 ? 1 : -1;
    int dz = cp_link->p()[2] - z*cpLength > cpLength/2 ? 1 : -1;
    // 接触点候補になりうる近傍のidxを求める
    std::vector<unsigned int> idxs;
    for (int kx=0; kx<2; kx++) {
      for (int ky=0; ky<2; ky++) {
        for (int kz=0; kz<2; kz++) {
          if (x+kx*dx >= cpAreaDim/2 || x+kx*dx < - cpAreaDim/2 ||
              y+ky*dy >= cpAreaDim/2 || y+ky*dy < - cpAreaDim/2 ||
              z+kz*dz >= cpAreaDim/2 || z+kz*dz < - cpAreaDim/2) continue;
          idxs.push_back(convertContactPointsIdx(x+kx*dx,y+ky*dy,z+kz*dz, cpAreaDim));
        }
      }
    }
    // 最近接接触点を選ぶ
    double minValue = 1e3;
    cnoid::Isometry3 selectedContact;
    for (unsigned int idx=0; idx<idxs.size(); idx++){
      for (int i=0; i<cps[idxs[idx]].size(); i++) {
        double value = (cps[idxs[idx]][i].translation() - cp_link->p()).norm();
        //          value -= std::min(this->weight_[3], this->weight_[4])*(cps[idxs[idx]][i].linear()*cnoid::Vector3::UnitZ()).dot(cp_link->R() * cnoid::Vector3::UnitZ());
        if (value < minValue) {
          minValue = value;
          selectedContact = cps[idxs[idx]][i];
          normal = normals[idxs[idx]][i];
          normalJacobian[0] = normalJXs[idxs[idx]][i];
          normalJacobian[1] = normalJYs[idxs[idx]][i];
          normalJacobian[2] = normalJZs[idxs[idx]][i];
        }
      }
    }
    if (minValue != 1e3) {
      localPos = selectedContact;
      cp_link->T() = selectedContact;
    } else { // 近傍にないときは前回の接触点を使う. 探索の安定化とglobal samplingのため
      cp_link->T() = localPos;
    }

  }

  void BodyContactConstraint::updateBounds() {
    if (this->A_contact_pos_link_ && !this->A_contact_pos_link_->isFreeJoint()) {
      std::cerr << "[BodyContactConstraint] A_contact_pos_link is not FreeJoint !" << std::endl;
    }

    if (this->B_contact_pos_link_ && !this->B_contact_pos_link_->isFreeJoint()) {
      std::cerr << "[BodyContactConstraint] B_contact_pos_link is not FreeJoint !" << std::endl;
    }

    if (this->A_contact_pos_link_ && (this->A_contactNormals_.size() == 0 ||
                                      this->A_contactNormalJacobianXs_.size() == 0 ||
                                      this->A_contactNormalJacobianYs_.size() == 0 ||
                                      this->A_contactNormalJacobianZs_.size() == 0)) {
      cnoid::TimeMeasure timer;
      if (this->debugLevel_ >= 1) timer.begin();
      calcNominals(this->A_contactPoints_, this->A_contactPointLength_, this->A_contactPointAreaDim_, this->A_contactNormals_, this->A_contactNormalJacobianXs_, this->A_contactNormalJacobianYs_, this->A_contactNormalJacobianZs_);
      if (this->debugLevel_ >= 1) std::cerr << "[BodyContactConstraint] " << (this->A_link_ ? this->A_link_->name() : std::string("world")) << " create nominals and jacobians" << timer.measure() << " [s]" << std::endl;
    }
    if (this->B_contact_pos_link_ && (this->B_contactNormals_.size() == 0 ||
                                      this->B_contactNormalJacobianXs_.size() == 0 ||
                                      this->B_contactNormalJacobianYs_.size() == 0 ||
                                      this->B_contactNormalJacobianZs_.size() == 0)) {
      cnoid::TimeMeasure timer;
      if (this->debugLevel_ >= 1) timer.begin();
      calcNominals(this->B_contactPoints_, this->B_contactPointLength_, this->B_contactPointAreaDim_, this->B_contactNormals_, this->B_contactNormalJacobianXs_, this->B_contactNormalJacobianYs_, this->B_contactNormalJacobianZs_);
      if (this->debugLevel_ >= 1) std::cerr << "[BodyContactConstraint] " << (this->B_link_ ? this->B_link_->name() : std::string("world")) << " create nominals and jacobians" << timer.measure() << " [s]" << std::endl;
    }
    // 探索された接触点位置をもとに、PositionConstraintのA_localpos_を更新する
    // 最も近い接触点候補を選ぶ
    if (this->A_contact_pos_link_) {
      cnoid::TimeMeasure timer;
      if (this->debugLevel_ >= 1) timer.begin();
      selectNearestPoint(this->A_contactPoints_, this->A_contactPointLength_, this->A_contactPointAreaDim_, this->A_contactNormals_, this->A_contactNormalJacobianXs_, this->A_contactNormalJacobianYs_, this->A_contactNormalJacobianZs_, A_contact_pos_link_, this->A_normal_, this->A_normalJacobian_, this->A_localpos_);
      if (this->debugLevel_ >= 1) std::cerr << "[BodyContactConstraint] " << (this->A_link_ ? this->A_link_->name() : std::string("world")) << " search nearest contact point " << timer.measure() << " [s]" << std::endl;
    }
    // 探索された接触点位置をもとに、PositionConstraintのB_localpos_を更新する
    // 最も近い接触点候補を選ぶ
    if (this->B_contact_pos_link_) {
      cnoid::TimeMeasure timer;
      if (this->debugLevel_ >= 1) timer.begin();
      selectNearestPoint(this->B_contactPoints_, this->B_contactPointLength_, this->B_contactPointAreaDim_, this->B_contactNormals_, this->B_contactNormalJacobianXs_, this->B_contactNormalJacobianYs_, this->B_contactNormalJacobianZs_, B_contact_pos_link_, this->B_normal_, this->B_normalJacobian_, this->B_localpos_);
      if (this->debugLevel_ >= 1) std::cerr << "[BodyContactConstraint] " << (this->B_link_ ? this->B_link_->name() : std::string("world")) << " search nearest contact point " << timer.measure() << " [s]" << std::endl;
    }
    if (this->debugLevel_ >= 2) {
      std::cerr << "selected A contact point" << std::endl;
      std::cerr << this->A_localpos_.translation().transpose() << std::endl;
      std::cerr << this->A_localpos_.rotation() << std::endl;
      std::cerr << "selected A contact normal : " << this->A_normal_.transpose() << std::endl;
      std::cerr << "selected A contact normal jacobian X : " << this->A_normalJacobian_[0].transpose() << std::endl;
      std::cerr << "selected A contact normal jacobian Y : " << this->A_normalJacobian_[1].transpose() << std::endl;
      std::cerr << "selected A contact normal jacobian Z : " << this->A_normalJacobian_[2].transpose() << std::endl;
      std::cerr << std::endl;
      std::cerr << "selected B contact point" << std::endl;
      std::cerr << this->B_localpos_.translation().transpose() << std::endl;
      std::cerr << this->B_localpos_.rotation() << std::endl;
      std::cerr << "selected B contact normal : " << this->B_normal_.transpose() << std::endl;
      std::cerr << "selected B contact normal jacobian X : " << this->B_normalJacobian_[0].transpose() << std::endl;
      std::cerr << "selected B contact normal jacobian Y : " << this->B_normalJacobian_[1].transpose() << std::endl;
      std::cerr << "selected B contact normal jacobian Z : " << this->B_normalJacobian_[2].transpose() << std::endl;
    }
    double prevDistance = this->distance();
    PositionConstraint::updateBounds(); // contactPointsの分解能によっては振動する可能性
    this->distanceDiff_ = std::abs(prevDistance - this->distance());
  }

  void BodyContactConstraint::updateJacobian (const std::vector<cnoid::LinkPtr>& joints) {

    // this->jacobian_を計算する
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!IKConstraint::isJointsSame(joints,this->jacobian_joints_)
       || this->A_link_ != this->jacobian_A_link_
       || this->B_link_ != this->jacobian_B_link_
       || this->eval_link_ != this->jacobian_eval_link_
       || this->A_contact_pos_link_ != this->A_jacobian_contact_pos_link_
       || this->B_contact_pos_link_ != this->B_jacobian_contact_pos_link_){
      this->jacobian_joints_ = joints;
      this->jacobian_A_link_ = this->A_link_;
      this->jacobian_B_link_ = this->B_link_;
      this->jacobian_eval_link_ = this->eval_link_;
      this->A_jacobian_contact_pos_link_ = this->A_contact_pos_link_;
      this->B_jacobian_contact_pos_link_ = this->B_contact_pos_link_;

      ik_constraint2::calc6DofJacobianShape(this->jacobian_joints_,//input
                                            this->jacobian_A_link_,//input
                                            this->jacobian_A_full_,
                                            this->jacobianColMap_,
                                            this->path_A_joints_
                                            );
      ik_constraint2::calc6DofJacobianShape(this->jacobian_joints_,//input
                                            this->jacobian_B_link_,//input
                                            this->jacobian_B_full_,
                                            this->jacobianColMap_,
                                            this->path_B_joints_
                                            );
      ik_constraint2::calc6DofJacobianShape(this->jacobian_joints_,//input
                                            this->jacobian_eval_link_,//input
                                            this->jacobian_eval_full_,
                                            this->jacobianColMap_,
                                            this->path_eval_joints_
                                            );
      if (this->jacobianColMap_.find(this->A_contact_pos_link_) != this->jacobianColMap_.end()) {
        int num_variables = 0;
        for(size_t i=0;i<joints.size();i++){
          num_variables += getJointDOF(joints[i]);
        }
        this->jacobian_A_full_.insert(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_A_full_.insert(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_A_full_.insert(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_A_full_.insert(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_A_full_.insert(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_A_full_.insert(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_A_full_.insert(2,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_A_full_.insert(2,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_A_full_.insert(2,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_A_full_.insert(3,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_A_full_.insert(3,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_A_full_.insert(3,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_A_full_.insert(4,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_A_full_.insert(4,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_A_full_.insert(4,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_A_full_.insert(5,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_A_full_.insert(5,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_A_full_.insert(5,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = 1;
        this->A_jacobian_contact_pos_= Eigen::SparseMatrix<double,Eigen::RowMajor>(1,num_variables);
        this->A_jacobian_contact_pos_.insert(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = 1;
        this->A_jacobian_contact_pos_.insert(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = 1;
        this->A_jacobian_contact_pos_.insert(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = 1;
        this->A_jacobian_ineq_contact_pos_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(2,num_variables);
        this->A_jacobian_ineq_contact_pos_.insert(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]) = 1;
        this->A_jacobian_ineq_contact_pos_.insert(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]) = 1;
      }

      if (this->jacobianColMap_.find(this->B_contact_pos_link_) != this->jacobianColMap_.end()) {
        int num_variables = 0;
        for(size_t i=0;i<joints.size();i++){
          num_variables += getJointDOF(joints[i]);
        }
        this->jacobian_B_full_.insert(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_B_full_.insert(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_B_full_.insert(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_B_full_.insert(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_B_full_.insert(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_B_full_.insert(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_B_full_.insert(2,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_B_full_.insert(2,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_B_full_.insert(2,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_B_full_.insert(3,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_B_full_.insert(3,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_B_full_.insert(3,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_B_full_.insert(4,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_B_full_.insert(4,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_B_full_.insert(4,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_B_full_.insert(5,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_B_full_.insert(5,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_B_full_.insert(5,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = 1;
        this->B_jacobian_contact_pos_= Eigen::SparseMatrix<double,Eigen::RowMajor>(1,num_variables);
        this->B_jacobian_contact_pos_.insert(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = 1;
        this->B_jacobian_contact_pos_.insert(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = 1;
        this->B_jacobian_contact_pos_.insert(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = 1;
        this->B_jacobian_ineq_contact_pos_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(2,num_variables);
        this->B_jacobian_ineq_contact_pos_.insert(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]) = 1;
        this->B_jacobian_ineq_contact_pos_.insert(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]) = 1;
      }
    }

    ik_constraint2::calc6DofJacobianCoef(this->jacobian_joints_,//input
                                         this->jacobian_A_link_,//input
                                         this->A_localpos_.translation(),//input
                                         this->jacobianColMap_,//input
                                         this->path_A_joints_,//input
                                         this->jacobian_A_full_
                                         );
    ik_constraint2::calc6DofJacobianCoef(this->jacobian_joints_,//input
                                         this->jacobian_B_link_,//input
                                         this->B_localpos_.translation(),//input
                                         this->jacobianColMap_,//input
                                         this->path_B_joints_,//input
                                         this->jacobian_B_full_
                                         );
    ik_constraint2::calc6DofJacobianCoef(this->jacobian_joints_,//input
                                         this->jacobian_eval_link_,//input
                                         cnoid::Vector3::Zero(),//input
                                         this->jacobianColMap_,//input
                                         this->path_eval_joints_,//input
                                         this->jacobian_eval_full_
                                         );

    if (this->jacobianColMap_.find(this->A_contact_pos_link_) != this->jacobianColMap_.end()) {
      this->jacobian_A_full_.coeffRef(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->A_link()->R()(0,0) : 0;
      this->jacobian_A_full_.coeffRef(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->A_link()->R()(0,1) : 0;
      this->jacobian_A_full_.coeffRef(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->A_link()->R()(0,2) : 0;
      this->jacobian_A_full_.coeffRef(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->A_link()->R()(1,0) : 0;
      this->jacobian_A_full_.coeffRef(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->A_link()->R()(1,1) : 0;
      this->jacobian_A_full_.coeffRef(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->A_link()->R()(1,2) : 0;
      this->jacobian_A_full_.coeffRef(2,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->A_link()->R()(2,0) : 0;
      this->jacobian_A_full_.coeffRef(2,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->A_link()->R()(2,1) : 0;
      this->jacobian_A_full_.coeffRef(2,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->A_link()->R()(2,2) : 0;
      cnoid::Matrix3d w;
      w << this->A_normalJacobian_[0][0], this->A_normalJacobian_[1][0], this->A_normalJacobian_[1][0],
           this->A_normalJacobian_[0][1], this->A_normalJacobian_[1][1], this->A_normalJacobian_[1][1],
           this->A_normalJacobian_[0][2], this->A_normalJacobian_[1][2], this->A_normalJacobian_[1][2];
      cnoid::Matrix3d Rw = this->A_link()->R() * w;
      this->jacobian_A_full_.coeffRef(3,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(0,0) : 0;
      this->jacobian_A_full_.coeffRef(3,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(0,1) : 0;
      this->jacobian_A_full_.coeffRef(3,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(0,2) : 0;
      this->jacobian_A_full_.coeffRef(4,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(1,0) : 0;
      this->jacobian_A_full_.coeffRef(4,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(1,1) : 0;
      this->jacobian_A_full_.coeffRef(4,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(1,2) : 0;
      this->jacobian_A_full_.coeffRef(5,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(2,0) : 0;
      this->jacobian_A_full_.coeffRef(5,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(2,1) : 0;
      this->jacobian_A_full_.coeffRef(5,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(2,2) : 0;
    }

    if (this->jacobianColMap_.find(this->B_contact_pos_link_) != this->jacobianColMap_.end()) {
      this->jacobian_B_full_.coeffRef(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->B_link()->R()(0,0) : 0;
      this->jacobian_B_full_.coeffRef(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->B_link()->R()(0,1) : 0;
      this->jacobian_B_full_.coeffRef(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->B_link()->R()(0,2) : 0;
      this->jacobian_B_full_.coeffRef(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->B_link()->R()(1,0) : 0;
      this->jacobian_B_full_.coeffRef(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->B_link()->R()(1,1) : 0;
      this->jacobian_B_full_.coeffRef(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->B_link()->R()(1,2) : 0;
      this->jacobian_B_full_.coeffRef(2,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->B_link()->R()(2,0) : 0;
      this->jacobian_B_full_.coeffRef(2,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->B_link()->R()(2,1) : 0;
      this->jacobian_B_full_.coeffRef(2,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * this->B_link()->R()(2,2) : 0;
      cnoid::Matrix3d w;
      w << this->B_normalJacobian_[0][0], this->B_normalJacobian_[1][0], this->B_normalJacobian_[1][0],
           this->B_normalJacobian_[0][1], this->B_normalJacobian_[1][1], this->B_normalJacobian_[1][1],
           this->B_normalJacobian_[0][2], this->B_normalJacobian_[1][2], this->B_normalJacobian_[1][2];
      cnoid::Matrix3d Rw = this->B_link()->R() * w;
      this->jacobian_B_full_.coeffRef(3,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(0,0) : 0;
      this->jacobian_B_full_.coeffRef(3,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(0,1) : 0;
      this->jacobian_B_full_.coeffRef(3,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(0,2) : 0;
      this->jacobian_B_full_.coeffRef(4,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(1,0) : 0;
      this->jacobian_B_full_.coeffRef(4,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(1,1) : 0;
      this->jacobian_B_full_.coeffRef(4,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(1,2) : 0;
      this->jacobian_B_full_.coeffRef(5,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(2,0) : 0;
      this->jacobian_B_full_.coeffRef(5,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(2,1) : 0;
      this->jacobian_B_full_.coeffRef(5,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = (this->distanceDiff_ < this->contactSearchLimit_) ? this->contactWeight_ * Rw(2,2) : 0;
    }

    cnoid::Matrix3d eval_R = (this->eval_link_) ? this->eval_link_->R() * this->eval_localR_ : this->eval_localR_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> eval_R_sparse(3,3);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) eval_R_sparse.insert(i,j) = eval_R(i,j);
    this->jacobian_full_local_.resize(6, this->jacobian_A_full_.cols());

    this->jacobian_full_local_.topRows<3>() = eval_R_sparse.transpose() * this->jacobian_A_full_.topRows<3>();
    this->jacobian_full_local_.bottomRows<3>() = eval_R_sparse.transpose() * this->jacobian_A_full_.bottomRows<3>();
    this->jacobian_full_local_.topRows<3>() -= Eigen::SparseMatrix<double,Eigen::RowMajor>(eval_R_sparse.transpose() * this->jacobian_B_full_.topRows<3>());
    this->jacobian_full_local_.bottomRows<3>() -= Eigen::SparseMatrix<double,Eigen::RowMajor>(eval_R_sparse.transpose() * this->jacobian_B_full_.bottomRows<3>());
    this->jacobian_full_local_.topRows<3>() += IKConstraint::cross(this->current_error_eval_.head<3>()) * eval_R_sparse.transpose() * this->jacobian_eval_full_.bottomRows<3>();
    this->jacobian_full_local_.bottomRows<3>() += IKConstraint::cross(this->current_error_eval_.tail<3>()) * eval_R_sparse.transpose() * this->jacobian_eval_full_.bottomRows<3>();

    this->jacobian_.resize((this->weight_.array() > 0.0).count(),this->jacobian_full_local_.cols());
    for(size_t i=0, idx=0;i<6;i++){
      if(this->weight_[i]>0.0) {
        this->jacobian_.row(idx) = this->weight_[i] * this->jacobian_full_local_.row(i);
        idx++;
      }
    }

    // 接触点の接平面でのみ動く
    if (this->jacobianColMap_.find(this->A_contact_pos_link_) != this->jacobianColMap_.end()) {
      cnoid::Vector3 normal = this->A_normal_;
      this->A_jacobian_contact_pos_.coeffRef(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = normal[0];
      this->A_jacobian_contact_pos_.coeffRef(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = normal[1];
      this->A_jacobian_contact_pos_.coeffRef(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = normal[2];

      cnoid::Vector3 tangent_x = (normal==cnoid::Vector3::UnitY()) ? cnoid::Vector3::UnitZ() : cnoid::Vector3::UnitY().cross(normal);
      this->A_jacobian_ineq_contact_pos_.coeffRef(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = tangent_x[0];
      this->A_jacobian_ineq_contact_pos_.coeffRef(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = tangent_x[1];
      this->A_jacobian_ineq_contact_pos_.coeffRef(0,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = tangent_x[2];
      cnoid::Vector3 tangent_y = normal.cross(tangent_x);
      this->A_jacobian_ineq_contact_pos_.coeffRef(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+0) = tangent_y[0];
      this->A_jacobian_ineq_contact_pos_.coeffRef(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+1) = tangent_y[1];
      this->A_jacobian_ineq_contact_pos_.coeffRef(1,this->jacobianColMap_[this->A_jacobian_contact_pos_link_]+2) = tangent_y[2];

      this->jacobianIneq_.resize(3,this->jacobian_.cols());
      if (this->distanceDiff_ < this->contactSearchLimit_) {
        this->jacobianIneq_.row(0) = this->A_jacobian_contact_pos_.row(0);
        this->jacobianIneq_.row(1) = this->A_jacobian_ineq_contact_pos_.row(0);
        this->jacobianIneq_.row(2) = this->A_jacobian_ineq_contact_pos_.row(1);
      }
      if (this->minIneq_.size() != 3 ||
          this->maxIneq_.size() != 3) {
        // 接触点の接平面でのみ動く
        // 接触点の1イテレーションあたり探索領域の制限
        this->minIneq_.resize(3);
        this->minIneq_[0] = 0.0;
        this->minIneq_[1] = - this->contactSearchLimit_;
        this->minIneq_[2] = - this->contactSearchLimit_;
        this->maxIneq_.resize(3);
        this->maxIneq_[0] = 0.0;
        this->maxIneq_[1] = this->contactSearchLimit_;
        this->maxIneq_[2] = this->contactSearchLimit_;
      }
    }

    if (this->jacobianColMap_.find(this->B_contact_pos_link_) != this->jacobianColMap_.end()) {
      cnoid::Vector3 normal = this->B_normal_;
      this->B_jacobian_contact_pos_.coeffRef(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = normal[0];
      this->B_jacobian_contact_pos_.coeffRef(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = normal[1];
      this->B_jacobian_contact_pos_.coeffRef(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = normal[2];

      cnoid::Vector3 tangent_x = (normal==cnoid::Vector3::UnitY()) ? cnoid::Vector3::UnitZ() : cnoid::Vector3::UnitY().cross(normal);
      this->B_jacobian_ineq_contact_pos_.coeffRef(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = tangent_x[0];
      this->B_jacobian_ineq_contact_pos_.coeffRef(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = tangent_x[1];
      this->B_jacobian_ineq_contact_pos_.coeffRef(0,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = tangent_x[2];
      cnoid::Vector3 tangent_y = normal.cross(tangent_x);
      this->B_jacobian_ineq_contact_pos_.coeffRef(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+0) = tangent_y[0];
      this->B_jacobian_ineq_contact_pos_.coeffRef(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+1) = tangent_y[1];
      this->B_jacobian_ineq_contact_pos_.coeffRef(1,this->jacobianColMap_[this->B_jacobian_contact_pos_link_]+2) = tangent_y[2];

      this->jacobianIneq_.resize(3,this->jacobian_.cols());
      if (this->distanceDiff_ < this->contactSearchLimit_) {
        this->jacobianIneq_.row(0) = this->B_jacobian_contact_pos_.row(0);
        this->jacobianIneq_.row(1) = this->B_jacobian_ineq_contact_pos_.row(0);
        this->jacobianIneq_.row(2) = this->B_jacobian_ineq_contact_pos_.row(1);
      }
      if (this->minIneq_.size() != 3 ||
          this->maxIneq_.size() != 3) {
        // 接触点の接平面でのみ動く
        // 接触点の1イテレーションあたり探索領域の制限
        this->minIneq_.resize(3);
        this->minIneq_[0] = 0.0;
        this->minIneq_[1] = - this->contactSearchLimit_;
        this->minIneq_[2] = - this->contactSearchLimit_;
        this->maxIneq_.resize(3);
        this->maxIneq_[0] = 0.0;
        this->maxIneq_[1] = this->contactSearchLimit_;
        this->maxIneq_[2] = this->contactSearchLimit_;
      }
    }

    if ((this->jacobianColMap_.find(this->A_contact_pos_link_) == this->jacobianColMap_.end()) && (this->jacobianColMap_.find(this->B_contact_pos_link_) == this->jacobianColMap_.end())) {
      // this->jacobianIneq_のサイズだけそろえる
      this->jacobianIneq_.resize(0,this->jacobian_.cols());
    }

    if(this->debugLevel_>=1){
      std::cerr << "BodyContactConstraint" << std::endl;
      std::cerr << "jacobian" << std::endl;
      std::cerr << this->jacobian_ << std::endl;
      std::cerr << "jacobianIneq" << std::endl;
      std::cerr << this->jacobianIneq_ << std::endl;
    }

    return;
  }

  std::shared_ptr<ik_constraint2::IKConstraint> BodyContactConstraint::clone(const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const {
    std::shared_ptr<BodyContactConstraint> ret = std::make_shared<BodyContactConstraint>(*this);
    this->copy(ret, modelMap);
    return ret;
  }

  void BodyContactConstraint::copy(std::shared_ptr<BodyContactConstraint> ret, const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const {
    if(this->A_link_ && modelMap.find(this->A_link_->body()) != modelMap.end()) ret->A_link() = modelMap.find(this->A_link_->body())->second->link(this->A_link_->index());
    if(this->B_link_ && modelMap.find(this->B_link_->body()) != modelMap.end()) ret->B_link() = modelMap.find(this->B_link_->body())->second->link(this->B_link_->index());
    if(this->eval_link_ && modelMap.find(this->eval_link_->body()) != modelMap.end()) ret->eval_link() = modelMap.find(this->eval_link_->body())->second->link(this->eval_link_->index());
    if(this->A_contact_pos_link_ && modelMap.find(this->A_contact_pos_link_->body()) != modelMap.end()) ret->A_contact_pos_link() = modelMap.find(this->A_contact_pos_link_->body())->second->link(this->A_contact_pos_link_->index());
    if(this->B_contact_pos_link_ && modelMap.find(this->B_contact_pos_link_->body()) != modelMap.end()) ret->B_contact_pos_link() = modelMap.find(this->B_contact_pos_link_->body())->second->link(this->B_contact_pos_link_->index());
  }

}
