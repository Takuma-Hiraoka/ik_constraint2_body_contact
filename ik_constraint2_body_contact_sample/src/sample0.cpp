#include <choreonoid_viewer/choreonoid_viewer.h>
#include <cnoid/Body>
#include <cnoid/MeshGenerator>
#include <prioritized_inverse_kinematics_solver2/prioritized_inverse_kinematics_solver2.h>
#include <ik_constraint2/ik_constraint2.h>
#include <ik_constraint2_body_contact/BodyContactConstraint.h>

namespace ik_constraint2_body_contact_sample{
  void sample0(){
    cnoid::MeshGenerator meshGenerator;
    cnoid::BodyPtr body = new cnoid::Body();
    cnoid::LinkPtr rootLink = new cnoid::Link();
    {
      cnoid::SgShapePtr shape = new cnoid::SgShape();
      shape->setMesh(meshGenerator.generateBox(cnoid::Vector3(1,1,0.1)));
      cnoid::SgMaterialPtr material = new cnoid::SgMaterial();
      material->setTransparency(0);
      shape->setMaterial(material);
      cnoid::SgPosTransformPtr posTransform = new cnoid::SgPosTransform();
      posTransform->translation() = cnoid::Vector3(0,0,-0.05);
      posTransform->addChild(shape);
      rootLink->addShapeNode(posTransform);
    }
    body->setRootLink(rootLink);
    cnoid::LinkPtr variable = new cnoid::Link();
    variable->setJointType(cnoid::Link::JointType::FreeJoint);

    // setup tasks
    std::vector<std::shared_ptr<ik_constraint2::IKConstraint> > constraints0;
    std::vector<std::shared_ptr<ik_constraint2::IKConstraint> > constraints1;
    {
      // task: rarm to target. never reach
      std::shared_ptr<ik_constraint2_body_contact::BodyContactConstraint> constraint = std::make_shared<ik_constraint2_body_contact::BodyContactConstraint>();
      constraint->A_link() = rootLink;
      constraint->B_link() = nullptr;
      constraint->B_localpos().translation() = cnoid::Vector3(0.3,0.3,0.0);
      constraint->contact_pos_link() = variable;
      constraint->debugLevel() = 0;
      double resolution=0.03;
      for (int i=0; i<1.0/resolution;i++) {
        for (int j=0; j<1.0/resolution;j++) {
          cnoid::Isometry3 contactPoint = cnoid::Isometry3::Identity();
          contactPoint.translation() = cnoid::Vector3(i*resolution - 0.5, j*resolution - 0.5, 0);
          constraint->contactPoints().push_back(contactPoint);
        }
      }
      constraint->precision() = resolution / 2;
      constraints1.push_back(constraint);
    }

    std::vector<std::shared_ptr<prioritized_qp_base::Task> > tasks;
    std::vector<cnoid::LinkPtr> variables{variable};
    std::vector<std::vector<std::shared_ptr<ik_constraint2::IKConstraint> > > constraints{constraints0, constraints1};
    std::shared_ptr<choreonoid_viewer::Viewer> viewer = std::make_shared<choreonoid_viewer::Viewer>();
    viewer->objects(body);

    viewer->drawObjects();
        // main loop
    for(int i=0;i<200;i++){
      prioritized_inverse_kinematics_solver2::IKParam param;
      param.debugLevel = 0; // debug
      param.maxIteration = 1;
      param.we = 1e2;
      bool solved = prioritized_inverse_kinematics_solver2::solveIKLoop(variables,
                                                                        constraints,
                                                                        tasks,
                                                                        param);

      if(i % 1 == 0){
        std::cerr << "loop: " << i << std::endl;
        std::vector<cnoid::SgNodePtr> markers;
        for(int j=0;j<constraints.size();j++){
          for(int k=0;k<constraints[j].size(); k++){
            const std::vector<cnoid::SgNodePtr>& marker = constraints[j][k]->getDrawOnObjects();
            std::copy(marker.begin(), marker.end(), std::back_inserter(markers));
          }
        }
        viewer->drawOn(markers);
        viewer->drawObjects();

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }

      if(solved) break;
    }

  }
}
