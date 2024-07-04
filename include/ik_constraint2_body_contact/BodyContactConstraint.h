#ifndef IK_CONSTRAINT2_BODYCONTACTCONSTRAINT_H
#define IK_CONSTRAINT2_BODYCONTACTCONSTRAINT_H

#include <ik_constraint2/ik_constraint2.h>

namespace ik_constraint2_body_contact{
  class BodyContactConstraint : public ik_constraint2::PositionConstraint {
  public:
    // A_linkの接触候補点の中から選んだA_localposとB_link中のB_localposを一致させる
    // 接触点の法線方向は考慮されるが、接平面の方向の自由度は制約しない. 法線方向の回転自由度の重みを0にすること.
    // contact_pos_link: 接触候補点探索用. variablesとして追加すること. このtranslationやrotationはA_linkのローカル座標系
    // contactPoints: 接触点候補
    const cnoid::LinkPtr& contact_pos_link() const { return contact_pos_link_;}
    cnoid::LinkPtr& contact_pos_link() { return contact_pos_link_;}
    const std::vector<cnoid::Isometry3>& contactPoints() const { return contactPoints_;}
    std::vector<cnoid::Isometry3>& contactPoints() { return contactPoints_;}
    // 探索された接触点位置をもとに、PositionConstraintのA_localpos_を更新する
    virtual void updateBounds() override;
    // 複製する. このとき、modelMapのkeyにあるロボットモデルに属するリンクは、valueに置き換える
    virtual std::shared_ptr<ik_constraint2::IKConstraint> clone(const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const override;
    void copy(std::shared_ptr<BodyContactConstraint> ret, const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const;

  private:
    cnoid::LinkPtr contact_pos_link_ = nullptr;
    std::vector<cnoid::Isometry3> contactPoints_;
  };
}

#endif
