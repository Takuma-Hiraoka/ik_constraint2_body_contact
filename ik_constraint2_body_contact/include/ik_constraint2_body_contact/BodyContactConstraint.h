#ifndef IK_CONSTRAINT2_BODYCONTACTCONSTRAINT_H
#define IK_CONSTRAINT2_BODYCONTACTCONSTRAINT_H

#include <ik_constraint2/ik_constraint2.h>

namespace ik_constraint2_body_contact{
  class BodyContactConstraint : public ik_constraint2::PositionConstraint {
  public:
    // A_linkの接触候補点の中から選んだA_localposとB_link中のB_localposを一致させる
    // 接触点の法線方向は考慮されるが、接平面の方向の自由度は制約しない. 法線方向の回転自由度の重みを0にすること.
    // precisionをcontactPointsの分解能と同程度にすること. そうしないとsatisfiedにならない.
    // contact_pos_link: 接触候補点探索用. variablesとして追加すること. JointはFreeJointにすること. このtranslationやrotationはA_linkのローカル座標系
    // contactPoints: 接触点候補
    const cnoid::LinkPtr& contact_pos_link() const { return contact_pos_link_;}
    cnoid::LinkPtr& contact_pos_link() { return contact_pos_link_;}
    const cnoid::BodyPtr& contact_pos_body() const { return contact_pos_body_;}
    cnoid::BodyPtr& contact_pos_body() { return contact_pos_body_;}
    const double& contactSearchLimit() const { return contactSearchLimit_;}
    double& contactSearchLimit() { return contactSearchLimit_;}
    const double& contactWeight() const { return contactWeight_;}
    double& contactWeight() { return contactWeight_;}
    const double& normalGradientDistance() const { return normalGradientDistance_;}
    double& normalGradientDistance() { return normalGradientDistance_;}
    // 探索された接触点位置をもとに、PositionConstraintのA_localpos_を更新する
    virtual void updateBounds() override;
    virtual void updateJacobian (const std::vector<cnoid::LinkPtr>& joints) override;

    unsigned int convertContactPointsIdx(int x, int y, int z);
    void setContactPoints(std::vector<cnoid::Isometry3> contactPoints, double contactPointLength=0.05, int contactPointAreaDim=16);
    // 複製する. このとき、modelMapのkeyにあるロボットモデルに属するリンクは、valueに置き換える
    virtual std::shared_ptr<ik_constraint2::IKConstraint> clone(const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const override;
    void copy(std::shared_ptr<BodyContactConstraint> ret, const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const;

  private:
    cnoid::LinkPtr contact_pos_link_ = nullptr;
    cnoid::BodyPtr contact_pos_body_ = nullptr; // clone用
    std::vector<std::vector<cnoid::Isometry3> > contactPoints_; // 計算速度向上のため、contactPointLength立方ごとに接触候補点をまとめておく.
    std::vector<std::vector<cnoid::Vector3> > contactNormals_; // contactPointとサイズが同じ. 接触点の法線方向を滑らかにしたもの. 角を通るヤコビアンを出すため.
    std::vector<std::vector<cnoid::Vector3> > contactNormalJacobianXs_; // contactPointとサイズが同じ. 接触点の近傍の法線をx偏微分したもの. 姿勢のヤコビアンに使う.
    std::vector<std::vector<cnoid::Vector3> > contactNormalJacobianYs_; // contactPointとサイズが同じ. 接触点の近傍の法線をy偏微分したもの. 姿勢のヤコビアンに使う.
    std::vector<std::vector<cnoid::Vector3> > contactNormalJacobianZs_; // contactPointとサイズが同じ. 接触点の近傍の法線をz偏微分したもの. 姿勢のヤコビアンに使う.
    double contactSearchLimit_ = 0.04; // 接触点探索は1イテレーションあたりこの距離以内で行う. 接触点探索を行うのは、接触点固定でのdistanceの変化量がcontactSearchLimit_より小さくなってから.
    double contactWeight_ = 1.0; // 関節角度に対すると接触点の探索重み. 小さいほうがより接触点を探索する. 接触点が変わるとIKが振動気味になるのは避けられないので、先に接触点を動かして最近傍点あたりに寄せてから関節角度で近づいていくとよい.
    cnoid::Vector3 normal_ = cnoid::Vector3::UnitX();
    std::vector<cnoid::Vector3> normalJacobian_ = std::vector<cnoid::Vector3>{cnoid::Vector3::Zero(), cnoid::Vector3::Zero(), cnoid::Vector3::Zero()};
    double normalGradientDistance_ = 0.05; // contactNormalsを計算する際に、normalGradientDistanceの範囲内のcontactPointsから計算する.
    double contactPointLength_ = 0.05; // contactPointAreaLength
    int  contactPointAreaDim_ = 16; // contactPointsは全て原点中心の一辺がcontactPointAreaLength*contactPointAreaDimの立方体の中に収まっている必要がある.
    double distanceDiff_ = 0.0;

    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_contact_pos_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_ineq_contact_pos_;
    cnoid::LinkPtr jacobian_contact_pos_link_ = nullptr;// 前回のjacobian計算時のcontact_pos_link
  };
}

#endif
