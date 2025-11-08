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
      constraint->B_localpos().translation() = cnoid::Vector3(0.47,-0.11,-0.1);
      constraint->B_localpos().linear() = cnoid::Matrix3(cnoid::AngleAxis(M_PI,cnoid::Vector3(0,1,0)));
      constraint->A_contact_pos_link() = variable;
      constraint->debugLevel() = 0;
      constraint->weight()[3] = 0.1;
      constraint->weight()[4] = 0.1;
      constraint->weight()[5] = 0;
      std::vector<cnoid::Isometry3> contactPoints;
      double resolution=0.02;
      for (int i=0; i<1.0/resolution;i++) {
        for (int j=0; j<1.0/resolution;j++) {
          cnoid::Isometry3 contactPoint = cnoid::Isometry3::Identity();
          contactPoint.translation() = cnoid::Vector3(i*resolution - 0.5+0.01, j*resolution - 0.5+0.01, 0);
          contactPoints.push_back(contactPoint);
        }
      }
      for (int i=0; i<1.0/resolution;i++) {
        for (int j=0; j<0.1/resolution;j++) {
          cnoid::Isometry3 contactPoint = cnoid::Isometry3::Identity();
          contactPoint.translation() = cnoid::Vector3(0.5, i*resolution - 0.5+0.01, -j*resolution-0.01);
          contactPoint.linear() = cnoid::Matrix3(cnoid::AngleAxis(M_PI / 2,cnoid::Vector3(0,1,0)));
          contactPoints.push_back(contactPoint);
        }
      }
      for (int i=0; i<1.0/resolution;i++) {
        for (int j=0; j<0.1/resolution;j++) {
          cnoid::Isometry3 contactPoint = cnoid::Isometry3::Identity();
          contactPoint.translation() = cnoid::Vector3(-0.5, i*resolution - 0.5+0.01, -j*resolution-0.01);
          contactPoint.linear() = cnoid::Matrix3(cnoid::AngleAxis(-M_PI / 2,cnoid::Vector3(0,1,0)));
          contactPoints.push_back(contactPoint);
        }
      }
      for (int i=0; i<1.0/resolution;i++) {
        for (int j=0; j<0.1/resolution;j++) {
          cnoid::Isometry3 contactPoint = cnoid::Isometry3::Identity();
          contactPoint.translation() = cnoid::Vector3(i*resolution - 0.5+0.01, 0.5, -j*resolution-0.01);
          contactPoint.linear() = cnoid::Matrix3(cnoid::AngleAxis(-M_PI / 2,cnoid::Vector3(1,0,0)));
          contactPoints.push_back(contactPoint);
        }
      }
      for (int i=0; i<1.0/resolution;i++) {
        for (int j=0; j<0.1/resolution;j++) {
          cnoid::Isometry3 contactPoint = cnoid::Isometry3::Identity();
          contactPoint.translation() = cnoid::Vector3(i*resolution - 0.5+0.01, -0.5, -j*resolution-0.01);
          contactPoint.linear() = cnoid::Matrix3(cnoid::AngleAxis(M_PI / 2,cnoid::Vector3(1,0,0)));
          contactPoints.push_back(contactPoint);
        }
      }
      for (int i=0; i<1.0/resolution;i++) {
        for (int j=0; j<1.0/resolution;j++) {
          cnoid::Isometry3 contactPoint = cnoid::Isometry3::Identity();
          contactPoint.translation() = cnoid::Vector3(i*resolution - 0.5+0.01, j*resolution - 0.5+0.01, -0.1);
          contactPoint.linear() = cnoid::Matrix3(cnoid::AngleAxis(M_PI,cnoid::Vector3(0,1,0)));
          contactPoints.push_back(contactPoint);
        }
      }
      constraint->setContactPointsA(contactPoints, 0.05, 28);
      constraint->precision() = resolution;
      constraints1.push_back(constraint);
    }

    std::vector<std::shared_ptr<prioritized_qp_base::Task> > tasks;
    std::vector<cnoid::LinkPtr> variables{variable};
    std::vector<std::vector<std::shared_ptr<ik_constraint2::IKConstraint> > > constraints{constraints0, constraints1};
    std::shared_ptr<choreonoid_viewer::Viewer> viewer = std::make_shared<choreonoid_viewer::Viewer>();
    viewer->objects(body);

    viewer->drawObjects();

    prioritized_inverse_kinematics_solver2::IKParam param;
    param.debugLevel = 3;
    param.maxIteration = 200;
    param.satisfiedConvergeLevel = 2; // 姿勢によって角から離れようとしない問題があるので収束しないが、こうしないと動き始めでBodyContactが動くより前にconvergeThreを満たしてしまう.
    param.we = 1e2;
    param.viewer = viewer;
    param.viewMilliseconds = -1;

    prioritized_inverse_kinematics_solver2::solveIKLoop(variables,
                                                        constraints,
                                                        tasks,
                                                        param);
  }
}
