#include <ik_constraint2_body_contact/BodyContactConstraint.h>
#include <ik_constraint2/Jacobian.h>

namespace ik_constraint2_body_contact{
  size_t getJointDOF(const cnoid::LinkPtr& joint) {
    if(joint->isRevoluteJoint() || joint->isPrismaticJoint()) return 1;
    else if(joint->isFreeJoint()) return 6;
    else return 0;
  }

  void BodyContactConstraint::updateBounds() {
    if (!this->contact_pos_link_ || !this->contact_pos_link_->isFreeJoint()) {
      std::cerr << "[BodyContactConstraint] contact_pos_link is not FreeJoint !" << std::endl;
    }

    if (this->contactNormals_.size() == 0) {
      double normalAngle = M_PI * 2 / 3; // 薄い板の場合の裏側は除く
      for (int i=0; i<contactPoints_.size(); i++) {
        cnoid::Vector3 normalSum = cnoid::Vector3::Zero();
        int weightSum = 0;
        for (int j=0; j<contactPoints_.size(); j++) {
          double dist = (contactPoints_[i].translation() - contactPoints_[j].translation()).norm();
          if (dist < this->normalGradientDistance_ &&
              (normalAngle >= std::acos(std::min(1.0,(std::max(-1.0,(contactPoints_[i].linear()*cnoid::Vector3::UnitZ()).dot(contactPoints_[j].linear() * cnoid::Vector3::UnitZ()))))))) {
            normalSum += contactPoints_[j].linear()*cnoid::Vector3::UnitZ();
            weightSum++;
          }
        }
        this->contactNormals_.push_back((normalSum / weightSum).normalized());
      }
    }
    // 探索された接触点位置をもとに、PositionConstraintのA_localpos_を更新する
    // 最も近い接触点候補を選ぶ
    double minValue = 1e3;
    cnoid::Isometry3 selectedContact;
    for (int i=0; i<this->contactPoints_.size(); i++) {
      double value = (this->contactPoints_[i].translation() - this->contact_pos_link_->p()).norm();
      if (value < minValue) {
        minValue = value;
        selectedContact = this->contactPoints_[i];
        this->normal_ = this->contactNormals_[i];
      }
    }
    if (minValue != 1e3) {
      this->A_localpos_ = selectedContact;
      this->contact_pos_link_->T() = selectedContact;
    }
    if (this->debugLevel_ >= 2) {
      std::cerr << "selected contact point" << std::endl;
      std::cerr << this->A_localpos_.translation().transpose() << std::endl;
      std::cerr << this->A_localpos_.rotation() << std::endl;
      std::cerr << "selected contact normal" << std::endl;
      std::cerr << this->normal_.transpose() << std::endl;
    }
    PositionConstraint::updateBounds(); // contactPointsの分解能によっては振動する可能性
    // 接触点の接平面でのみ動く
    if (this->contact_pos_link_) {
      this->eq_.conservativeResize(this->eq_.size()+1);
      this->eq_[this->eq_.size()-1] = 0.0;
      if (this->minIneq_.size() != 2 ||
          this->maxIneq_.size() != 2) {
        // 接触点の1イテレーションあたり探索領域の制限
        this->minIneq_.resize(2);
        this->minIneq_[0] = - this->contactSearchLimit_;
        this->minIneq_[1] = - this->contactSearchLimit_;
        this->maxIneq_.resize(2);
        this->maxIneq_[0] = this->contactSearchLimit_;
        this->maxIneq_[1] = this->contactSearchLimit_;
      }
    }
  }

  void BodyContactConstraint::updateJacobian (const std::vector<cnoid::LinkPtr>& joints) {

    // this->jacobian_を計算する
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!IKConstraint::isJointsSame(joints,this->jacobian_joints_)
       || this->A_link_ != this->jacobian_A_link_
       || this->B_link_ != this->jacobian_B_link_
       || this->eval_link_ != this->jacobian_eval_link_
       || this->contact_pos_link_ != this->jacobian_contact_pos_link_){
      this->jacobian_joints_ = joints;
      this->jacobian_A_link_ = this->A_link_;
      this->jacobian_B_link_ = this->B_link_;
      this->jacobian_eval_link_ = this->eval_link_;
      this->jacobian_contact_pos_link_ = this->contact_pos_link_;

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
      if (this->jacobianColMap_.find(this->contact_pos_link_) != this->jacobianColMap_.end()) {
        int num_variables = 0;
        for(size_t i=0;i<joints.size();i++){
          num_variables += getJointDOF(joints[i]);
        }
        this->jacobian_A_full_.insert(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_A_full_.insert(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_A_full_.insert(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_A_full_.insert(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_A_full_.insert(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_A_full_.insert(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_A_full_.insert(2,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_A_full_.insert(2,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_A_full_.insert(2,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_contact_pos_= Eigen::SparseMatrix<double,Eigen::RowMajor>(1,num_variables);
        this->jacobian_contact_pos_.insert(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = 1;
        this->jacobian_contact_pos_.insert(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = 1;
        this->jacobian_contact_pos_.insert(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = 1;
        this->jacobian_ineq_contact_pos_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(2,num_variables);
        this->jacobian_ineq_contact_pos_.insert(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]) = 1;
        this->jacobian_ineq_contact_pos_.insert(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]) = 1;
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

    if (this->jacobianColMap_.find(this->contact_pos_link_) != this->jacobianColMap_.end()) {
      this->jacobian_A_full_.coeffRef(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = this->A_link()->R()(0,0);
      this->jacobian_A_full_.coeffRef(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = this->A_link()->R()(0,1);
      this->jacobian_A_full_.coeffRef(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = this->A_link()->R()(0,2);
      this->jacobian_A_full_.coeffRef(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = this->A_link()->R()(1,0);
      this->jacobian_A_full_.coeffRef(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = this->A_link()->R()(1,1);
      this->jacobian_A_full_.coeffRef(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = this->A_link()->R()(1,2);
      this->jacobian_A_full_.coeffRef(2,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = this->A_link()->R()(2,0);
      this->jacobian_A_full_.coeffRef(2,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = this->A_link()->R()(2,1);
      this->jacobian_A_full_.coeffRef(2,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = this->A_link()->R()(2,2);
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

    // 接触点の接平面でのみ動く
    int normaleq = 0;
    if (this->jacobianColMap_.find(this->contact_pos_link_) != this->jacobianColMap_.end()) normaleq = 1;
    this->jacobian_.resize((this->weight_.array() > 0.0).count()+normaleq,this->jacobian_full_local_.cols());
    for(size_t i=0, idx=0;i<6;i++){
      if(this->weight_[i]>0.0) {
        this->jacobian_.row(idx) = this->weight_[i] * this->jacobian_full_local_.row(i);
        idx++;
      }
    }

    if (this->jacobianColMap_.find(this->contact_pos_link_) != this->jacobianColMap_.end()) {
      cnoid::Vector3 normal = this->normal_;
      this->jacobian_contact_pos_.coeffRef(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = normal[0];
      this->jacobian_contact_pos_.coeffRef(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = normal[1];
      this->jacobian_contact_pos_.coeffRef(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = normal[2];
      this->jacobian_.row((this->weight_.array() > 0.0).count()) = this->jacobian_contact_pos_.row(0);

      cnoid::Vector3 tangent_x = (normal==cnoid::Vector3::UnitY()) ? cnoid::Vector3::UnitZ() : cnoid::Vector3::UnitY().cross(normal);
      this->jacobian_ineq_contact_pos_.coeffRef(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = tangent_x[0];
      this->jacobian_ineq_contact_pos_.coeffRef(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = tangent_x[1];
      this->jacobian_ineq_contact_pos_.coeffRef(0,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = tangent_x[2];
      cnoid::Vector3 tangent_y = normal.cross(tangent_x);
      this->jacobian_ineq_contact_pos_.coeffRef(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]+0) = tangent_y[0];
      this->jacobian_ineq_contact_pos_.coeffRef(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]+1) = tangent_y[1];
      this->jacobian_ineq_contact_pos_.coeffRef(1,this->jacobianColMap_[this->jacobian_contact_pos_link_]+2) = tangent_y[2];

      this->jacobianIneq_.resize(2,this->jacobian_.cols());
      this->jacobianIneq_.row(0) = this->jacobian_ineq_contact_pos_.row(0);
      this->jacobianIneq_.row(1) = this->jacobian_ineq_contact_pos_.row(1);
    } else {
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
  }

}
