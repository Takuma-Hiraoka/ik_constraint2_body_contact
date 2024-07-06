#include <choreonoid_viewer/choreonoid_viewer.h>
#include <cnoid/Body>
#include <cnoid/BodyLoader>
#include <cnoid/SceneMarkers>
#include <cnoid/YAMLReader>
#include <iostream>
#include <ros/package.h>
#include <prioritized_inverse_kinematics_solver2/prioritized_inverse_kinematics_solver2.h>
#include <ik_constraint2/ik_constraint2.h>
#include <ik_constraint2_body_contact/BodyContactConstraint.h>

namespace ik_constraint2_body_contact_sample{
  void sample1(){
    cnoid::BodyLoader bodyLoader;
    cnoid::BodyPtr robot = bodyLoader.load(ros::package::getPath("choreonoid") + "/share/model/SR1/SR1.body");
    robot->rootLink()->p() = cnoid::Vector3(0,0,0.6);
    robot->rootLink()->v().setZero();
    robot->rootLink()->R() = cnoid::Matrix3::Identity();
    robot->rootLink()->w().setZero();
    std::vector<double> reset_manip_pose{
      0.0, -0.349066, 0.0, 0.820305, -0.471239, 0.0,// rleg
        0.523599, 0.0, 0.0, -1.74533, 0.15708, -0.113446, 0.637045,// rarm
        0.0, -0.349066, 0.0, 0.820305, -0.471239, 0.0,// lleg
        0.523599, 0.0, 0.0, -1.74533, -0.15708, -0.113446, -0.637045,// larm
        0.0, 0.0, 0.0};
    for(int j=0; j < robot->numJoints(); ++j){
      robot->joint(j)->q() = reset_manip_pose[j];
    }
    robot->calcForwardKinematics();
    robot->calcCenterOfMass();

    std::string contactFileName = ros::package::getPath("ik_constraint2_body_contact_sample") + "/config/sample_config.yaml";
    cnoid::YAMLReader reader;
    cnoid::MappingPtr node;
    std::string prevLinkName = "";
    std::vector<cnoid::Isometry3> contactPoints;
    try {
      node = reader.loadDocument(contactFileName)->toMapping();
    } catch(const cnoid::ValueNode::Exception& ex) {
      std::cerr << ex.message()  << std::endl;
    }
    if(node){
      cnoid::Listing* tactileSensorList = node->findListing("tactile_sensor");
      if (!tactileSensorList->isValid()) {
        std::cerr << "tactile_sensor list is not valid" << std::endl;
      }else{
        for (int i=0; i< tactileSensorList->size(); i++) {
          cnoid::Mapping* info = tactileSensorList->at(i)->toMapping();
          std::string linkName;
          // linkname
          info->extract("link", linkName);
          if (linkName != "RARM_WRIST_R") continue;
          cnoid::Isometry3 sensor;
          // translation
          cnoid::ValueNodePtr translation_ = info->extract("translation");
          if(translation_){
            cnoid::ListingPtr translationTmp = translation_->toListing();
            if(translationTmp->size()==3){
              sensor.translation() = cnoid::Vector3(translationTmp->at(0)->toDouble(), translationTmp->at(1)->toDouble(), translationTmp->at(2)->toDouble());
            }
          }
          // rotation
          cnoid::ValueNodePtr rotation_ = info->extract("rotation");
          if(rotation_){
            cnoid::ListingPtr rotationTmp = rotation_->toListing();
            if(rotationTmp->size() == 4){
              sensor.linear() = cnoid::AngleAxisd(rotationTmp->at(3)->toDouble(),
                                                  cnoid::Vector3{rotationTmp->at(0)->toDouble(), rotationTmp->at(1)->toDouble(), rotationTmp->at(2)->toDouble()}).toRotationMatrix();
            }
          }
          contactPoints.push_back(sensor);
        }
      }
    }


    // setup tasks
    std::vector<std::shared_ptr<ik_constraint2::IKConstraint> > constraints0;
    // joint limit
    for(int i=0;i<robot->numJoints();i++){
      std::shared_ptr<ik_constraint2::JointLimitConstraint> constraint = std::make_shared<ik_constraint2::JointLimitConstraint>();
      constraint->joint() = robot->joint(i);
      constraints0.push_back(constraint);
    }
    std::vector<std::shared_ptr<ik_constraint2::IKConstraint> > constraints1;
    {
      // task: rleg to target
      std::shared_ptr<ik_constraint2::PositionConstraint> constraint = std::make_shared<ik_constraint2::PositionConstraint>();
      constraint->A_link() = robot->link("RLEG_ANKLE_R");
      constraint->A_localpos().translation() = cnoid::Vector3(0.0,0.0,-0.04);
      constraint->B_link() = nullptr;
      constraint->B_localpos().translation() = cnoid::Vector3(0.0,-0.2,-0.0);
      constraints1.push_back(constraint);
    }
    {
      // task: lleg to target
      std::shared_ptr<ik_constraint2::PositionConstraint> constraint = std::make_shared<ik_constraint2::PositionConstraint>();
      constraint->A_link() = robot->link("LLEG_ANKLE_R");
      constraint->A_localpos().translation() = cnoid::Vector3(0.0,0.0,-0.04);
      constraint->B_link() = nullptr;
      constraint->B_localpos().translation() = cnoid::Vector3(0.0,0.2,0.0);
      constraints1.push_back(constraint);
    }
    {
      // task: COM to target
      std::shared_ptr<ik_constraint2::COMConstraint> constraint = std::make_shared<ik_constraint2::COMConstraint>();
      constraint->A_robot() = robot;
      constraint->B_localp() = cnoid::Vector3(0.0,0.0,0.6);
      constraints1.push_back(constraint);
    }
    std::vector<std::shared_ptr<ik_constraint2::IKConstraint> > constraints2;
    cnoid::LinkPtr variable = new cnoid::Link();
    variable->setJointType(cnoid::Link::JointType::FreeJoint);
    {
      // task: rarm to target. never reach
      std::shared_ptr<ik_constraint2_body_contact::BodyContactConstraint> constraint = std::make_shared<ik_constraint2_body_contact::BodyContactConstraint>();
      constraint->A_link() = robot->link("RARM_WRIST_R");
      constraint->A_localpos().translation() = cnoid::Vector3(0.025,0.0,0.-0.01);
      constraint->B_link() = nullptr;
      constraint->B_localpos().translation() = cnoid::Vector3(0.85,-0.2,0.8);
      constraint->B_localpos().linear() = cnoid::Matrix3(cnoid::AngleAxis(M_PI/2,cnoid::Vector3(0,1,0)));
      constraint->eval_localR() = constraint->B_localpos().linear();
      constraint->contact_pos_link() = variable;
      constraint->contact_pos_link()->T() = constraint->A_localpos();
      constraint->setContactPoints(contactPoints, 0.05, 10);
      constraint->contactSearchLimit() = 0.02;
      constraint->precision() = 0.02;
      constraint->contactWeight() = 0.5;
      constraint->weight()[3] = 0.1;
      constraint->weight()[4] = 0.1;
      constraint->weight()[5] = 0;
      constraint->debugLevel() = 0;

      constraints2.push_back(constraint);
    }

    std::vector<std::shared_ptr<prioritized_qp_base::Task> > tasks;
    std::vector<cnoid::LinkPtr> variables;
    variables.push_back(robot->rootLink());
    for(size_t i=0;i<robot->numJoints();i++){
      variables.push_back(robot->joint(i));
    }
    variables.push_back(variable);
    std::vector<std::vector<std::shared_ptr<ik_constraint2::IKConstraint> > > constraints{constraints0,constraints1,constraints2};
    // setup viewer
    std::shared_ptr<choreonoid_viewer::Viewer> viewer = std::make_shared<choreonoid_viewer::Viewer>();
    viewer->objects(robot);

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
        getchar();
      }

      if(solved) break;
    }

    for(size_t i=0;i<constraints.size();i++){
      for(size_t j=0;j<constraints[i].size();j++){
        constraints[i][j]->debugLevel() = 0;//not debug
        constraints[i][j]->updateBounds();
        if(constraints[i][j]->isSatisfied()) std::cerr << "constraint " << i << " " << j << ": satisfied"<< std::endl;
        else std::cerr << "constraint " << i << " " << j << ": NOT satidfied"<< std::endl;
      }
    }

  }
}
