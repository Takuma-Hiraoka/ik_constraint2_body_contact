#include <ik_constraint2_body_contact/BodyContactConstraint.h>

namespace ik_constraint2_body_contact{
  void BodyContactConstraint::updateBounds() {
    // 探索された接触点位置をもとに、PositionConstraintのA_localpos_を更新する
    // 最も近い接触点候補を選ぶ
    double maxValue = -1e3;
    double weight = 0.1;
    cnoid::Isometry3 selectedContact;
    for (int i=0; i<this->contactPoints_.size(); i++) {
      double value=weight*(this->contactPoints_[i].linear()*cnoid::Vector3::UnitZ()).dot(this->contact_pos_link_->R() * cnoid::Vector3::UnitZ());
      value -= (this->contactPoints_[i].translation() - this->contact_pos_link_->p()).norm();
      if (value > maxValue) {
        maxValue = value;
        selectedContact = this->contactPoints_[i];
      }
    }
    this->A_localpos() = selectedContact;
    PositionConstraint::updateBounds();
  }
  std::shared_ptr<ik_constraint2::IKConstraint> BodyContactConstraint::clone(const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const {
    std::shared_ptr<BodyContactConstraint> ret = std::make_shared<BodyContactConstraint>(*this);
    this->copy(ret, modelMap);
    return ret;
  }

  void BodyContactConstraint::copy(std::shared_ptr<BodyContactConstraint> ret, const std::map<cnoid::BodyPtr, cnoid::BodyPtr>& modelMap) const {
    if(this->A_link() && modelMap.find(this->A_link()->body()) != modelMap.end()) ret->A_link() = modelMap.find(this->A_link()->body())->second->link(this->A_link()->index());
    if(this->B_link() && modelMap.find(this->B_link()->body()) != modelMap.end()) ret->B_link() = modelMap.find(this->B_link()->body())->second->link(this->B_link()->index());
    if(this->eval_link() && modelMap.find(this->eval_link()->body()) != modelMap.end()) ret->eval_link() = modelMap.find(this->eval_link()->body())->second->link(this->eval_link()->index());
  }

}
