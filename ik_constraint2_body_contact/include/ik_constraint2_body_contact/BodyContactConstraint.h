#ifndef IK_CONSTRAINT2_BODYCONTACTCONSTRAINT_H
#define IK_CONSTRAINT2_BODYCONTACTCONSTRAINT_H

#include <ik_constraint2/ik_constraint2.h>

namespace ik_constraint2_body_contact{
  class BodyContactConstraint : public ik_constraint2::PositionConstraint {
  public:
    // A_linkの接触候補点の中から選んだA_localposとB_linkの接触候補店から選んだB_localposを一致させる
    // 接触点の法線方向は考慮されるが、接平面の方向の自由度は制約しない. 法線方向の回転自由度の重みを0にすること.
    // precisionをcontactPointsの分解能と同程度にすること. そうしないとsatisfiedにならない.
    // *_contact_pos_link: 接触候補点探索用. variablesとして追加すること. JointはFreeJointにすること. このtranslationやrotationはA_linkのローカル座標系
    // *_contact_pos_body: clone用. *_contact_pos_linkをリンクとして持つ. clone時に必ずmodelMapに含めるようにすること. contact_pos_linkと探索変数に加えるvariableを一致させるため
    // *_contactPoints: 接触点候補
    const cnoid::LinkPtr& A_contact_pos_link() const { return A_contact_pos_link_;}
    cnoid::LinkPtr& A_contact_pos_link() { return A_contact_pos_link_;}
    const cnoid::BodyPtr& A_contact_pos_body() const { return A_contact_pos_body_;}
    cnoid::BodyPtr& A_contact_pos_body() { return A_contact_pos_body_;}
    const cnoid::LinkPtr& B_contact_pos_link() const { return B_contact_pos_link_;}
    cnoid::LinkPtr& B_contact_pos_link() { return B_contact_pos_link_;}
    const cnoid::BodyPtr& B_contact_pos_body() const { return B_contact_pos_body_;}
    cnoid::BodyPtr& B_contact_pos_body() { return B_contact_pos_body_;}
    const std::vector<std::vector<cnoid::Isometry3> >& A_contactPoints() const { return A_contactPoints_;}
    std::vector<std::vector<cnoid::Isometry3> >& A_contactPoints() { return A_contactPoints_;}
    const std::vector<std::vector<cnoid::Isometry3> >& B_contactPoints() const { return B_contactPoints_;}
    std::vector<std::vector<cnoid::Isometry3> >& B_contactPoints() { return B_contactPoints_;}
    const std::vector<std::vector<cnoid::Vector3> >& A_contactNormals() const { return A_contactNormals_;}
    std::vector<std::vector<cnoid::Vector3> >& A_contactNormals() { return A_contactNormals_;}
    const std::vector<std::vector<cnoid::Vector3> >& A_contactNormalJacobianXs() const { return A_contactNormalJacobianXs_;}
    std::vector<std::vector<cnoid::Vector3> >& A_contactNormalJacobianXs() { return A_contactNormalJacobianXs_;}
    const std::vector<std::vector<cnoid::Vector3> >& A_contactNormalJacobianYs() const { return A_contactNormalJacobianYs_;}
    std::vector<std::vector<cnoid::Vector3> >& A_contactNormalJacobianYs() { return A_contactNormalJacobianYs_;}
    const std::vector<std::vector<cnoid::Vector3> >& A_contactNormalJacobianZs() const { return A_contactNormalJacobianZs_;}
    std::vector<std::vector<cnoid::Vector3> >& A_contactNormalJacobianZs() { return A_contactNormalJacobianZs_;}
    const std::vector<std::vector<cnoid::Vector3> >& B_contactNormals() const { return B_contactNormals_;}
    std::vector<std::vector<cnoid::Vector3> >& B_contactNormals() { return B_contactNormals_;}
    const std::vector<std::vector<cnoid::Vector3> >& B_contactNormalJacobianXs() const { return B_contactNormalJacobianXs_;}
    std::vector<std::vector<cnoid::Vector3> >& B_contactNormalJacobianXs() { return B_contactNormalJacobianXs_;}
    const std::vector<std::vector<cnoid::Vector3> >& B_contactNormalJacobianYs() const { return B_contactNormalJacobianYs_;}
    std::vector<std::vector<cnoid::Vector3> >& B_contactNormalJacobianYs() { return B_contactNormalJacobianYs_;}
    const std::vector<std::vector<cnoid::Vector3> >& B_contactNormalJacobianZs() const { return B_contactNormalJacobianZs_;}
    std::vector<std::vector<cnoid::Vector3> >& B_contactNormalJacobianZs() { return B_contactNormalJacobianZs_;}
    const double& contactSearchLimit() const { return contactSearchLimit_;}
    double& contactSearchLimit() { return contactSearchLimit_;}
    const double& contactWeight() const { return contactWeight_;}
    double& contactWeight() { return contactWeight_;}
    const double& normalGradientDistance() const { return normalGradientDistance_;}
    double& normalGradientDistance() { return normalGradientDistance_;}
    // 探索された接触点位置をもとに、PositionConstraintのA_localpos_を更新する
    virtual void updateBounds() override;
    virtual void updateJacobian (const std::vector<cnoid::LinkPtr>& joints) override;

    unsigned int convertContactPointsIdx(int x, int y, int z, int dim);
    void A_setContactPoints(std::vector<cnoid::Isometry3> contactPoints, double contactPointLength=0.05, int contactPointAreaDim=16);
    void B_setContactPoints(std::vector<cnoid::Isometry3> contactPoints, double contactPointLength=0.05, int contactPointAreaDim=16);
    void calcNominals(const std::vector<std::vector<cnoid::Isometry3> >& cps, // input
                      const double& cpLength,
                      const int& cpAreaDim,
                      std::vector<std::vector<cnoid::Vector3> >& normals, // output
                      std::vector<std::vector<cnoid::Vector3> >& normalJXs,
                      std::vector<std::vector<cnoid::Vector3> >& normalJYs,
                      std::vector<std::vector<cnoid::Vector3> >& normalJZs);
    void selectNearestPoint(const std::vector<std::vector<cnoid::Isometry3> >& cps, // input
                            const double& cpLength,
                            const double& cpAreaDim,
                            const std::vector<std::vector<cnoid::Vector3> >& normals,
                            const std::vector<std::vector<cnoid::Vector3> >& normalJXs,
                            const std::vector<std::vector<cnoid::Vector3> >& normalJYs,
                            const std::vector<std::vector<cnoid::Vector3> >& normalJZs,
                            const cnoid::LinkPtr cp_link, // output
                            cnoid::Vector3& normal,
                            std::vector<cnoid::Vector3>& normalJ,
                            cnoid::Isometry3& localPos);
    // 複製する. このとき、modelMapのkeyにあるロボットモデルに属するリンクは、valueに置き換える
    virtual std::shared_ptr<ik_constraint2::IKConstraint> clone(const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const override;
    void copy(std::shared_ptr<BodyContactConstraint> ret, const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const;

  private:
    cnoid::LinkPtr A_contact_pos_link_ = nullptr;
    cnoid::BodyPtr A_contact_pos_body_ = nullptr;
    std::vector<std::vector<cnoid::Isometry3> > A_contactPoints_; // 計算速度向上のため、A_contactPointLength立方ごとに接触候補点をまとめておく.
    std::vector<std::vector<cnoid::Vector3> > A_contactNormals_; // A_contactPointとサイズが同じ. 接触点の法線方向を滑らかにしたもの. 角を通るヤコビアンを出すため.
    std::vector<std::vector<cnoid::Vector3> > A_contactNormalJacobianXs_; // contactPointとサイズが同じ. 接触点の近傍の法線をx偏微分したもの. 姿勢のヤコビアンに使う.
    std::vector<std::vector<cnoid::Vector3> > A_contactNormalJacobianYs_; // contactPointとサイズが同じ. 接触点の近傍の法線をy偏微分したもの. 姿勢のヤコビアンに使う.
    std::vector<std::vector<cnoid::Vector3> > A_contactNormalJacobianZs_; // contactPointとサイズが同じ. 接触点の近傍の法線をz偏微分したもの. 姿勢のヤコビアンに使う.
    cnoid::Vector3 A_normal_ = cnoid::Vector3::UnitX();
    std::vector<cnoid::Vector3> A_normalJacobian_ = std::vector<cnoid::Vector3>{cnoid::Vector3::Zero(), cnoid::Vector3::Zero(), cnoid::Vector3::Zero()};
    double A_contactPointLength_ = 0.05; // contactPointAreaLength
    int  A_contactPointAreaDim_ = 16; // contactPointsは全て原点中心の一辺がcontactPointAreaLength*contactPointAreaDimの立方体の中に収まっている必要がある.

    cnoid::LinkPtr B_contact_pos_link_ = nullptr;
    cnoid::BodyPtr B_contact_pos_body_ = nullptr;
    std::vector<std::vector<cnoid::Isometry3> > B_contactPoints_; // 計算速度向上のため、B_contactPointLength立方ごとに接触候補点をまとめておく.
    std::vector<std::vector<cnoid::Vector3> > B_contactNormals_; // B_contactPointとサイズが同じ. 接触点の法線方向を滑らかにしたもの. 角を通るヤコビアンを出すため.
    std::vector<std::vector<cnoid::Vector3> > B_contactNormalJacobianXs_; // contactPointとサイズが同じ. 接触点の近傍の法線をx偏微分したもの. 姿勢のヤコビアンに使う.
    std::vector<std::vector<cnoid::Vector3> > B_contactNormalJacobianYs_; // contactPointとサイズが同じ. 接触点の近傍の法線をy偏微分したもの. 姿勢のヤコビアンに使う.
    std::vector<std::vector<cnoid::Vector3> > B_contactNormalJacobianZs_; // contactPointとサイズが同じ. 接触点の近傍の法線をz偏微分したもの. 姿勢のヤコビアンに使う.
    cnoid::Vector3 B_normal_ = cnoid::Vector3::UnitX();
    std::vector<cnoid::Vector3> B_normalJacobian_ = std::vector<cnoid::Vector3>{cnoid::Vector3::Zero(), cnoid::Vector3::Zero(), cnoid::Vector3::Zero()};
    double B_contactPointLength_ = 0.05; // contactPointAreaLength
    int  B_contactPointAreaDim_ = 16; // contactPointsは全て原点中心の一辺がcontactPointAreaLength*contactPointAreaDimの立方体の中に収まっている必要がある.

    double contactSearchLimit_ = 0.04; // 接触点探索は1イテレーションあたりこの距離以内で行う. 接触点探索を行うのは、接触点固定でのdistanceの変化量がcontactSearchLimit_より小さくなってから.
    double contactWeight_ = 1.0; // 関節角度に対すると接触点の探索重み. 小さいほうがより接触点を探索する. 接触点が変わるとIKが振動気味になるのは避けられないので、先に接触点を動かして最近傍点あたりに寄せてから関節角度で近づいていくとよい.
    double normalGradientDistance_ = 0.05; // contactNormalsを計算する際に、normalGradientDistanceの範囲内のcontactPointsから計算する.
    double distanceDiff_ = 0.0;

    Eigen::SparseMatrix<double,Eigen::RowMajor> A_jacobian_contact_pos_; // 接平面固定
    Eigen::SparseMatrix<double,Eigen::RowMajor> A_jacobian_ineq_contact_pos_; // 移動量制限
    cnoid::LinkPtr A_jacobian_contact_pos_link_ = nullptr;// 前回のjacobian計算時のA_contact_pos_link

    Eigen::SparseMatrix<double,Eigen::RowMajor> B_jacobian_contact_pos_; // 接平面固定
    Eigen::SparseMatrix<double,Eigen::RowMajor> B_jacobian_ineq_contact_pos_; // 移動量制限
    cnoid::LinkPtr B_jacobian_contact_pos_link_ = nullptr;// 前回のjacobian計算時のB_contact_pos_link
};
}

#endif
