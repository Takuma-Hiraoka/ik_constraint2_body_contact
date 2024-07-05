#include <cnoid/Plugin>
#include <cnoid/ItemManager>

#include <choreonoid_viewer/choreonoid_viewer.h>

namespace ik_constraint2_body_contact_sample{
  void sample0();
  class sample0Item : public choreonoid_viewer::ViewerBaseItem {
  public:
    static void initializeClass(cnoid::ExtensionManager* ext){ ext->itemManager().registerClass<sample0Item>("sample0Item"); }
  protected:
    virtual void main() override{ sample0(); return;}
  };
  typedef cnoid::ref_ptr<sample0Item> sample0ItemPtr;

  void sample1();
  class sample1Item : public choreonoid_viewer::ViewerBaseItem {
  public:
    static void initializeClass(cnoid::ExtensionManager* ext){ ext->itemManager().registerClass<sample1Item>("sample1Item"); }
  protected:
    virtual void main() override{ sample1(); return;}
  };
  typedef cnoid::ref_ptr<sample1Item> sample1ItemPtr;

  class IKConstraintBodyContactSamplePlugin : public cnoid::Plugin
  {
  public:
    IKConstraintBodyContactSamplePlugin() : Plugin("IKConstraintBodyContactSample")
    {
      require("Body");
    }
    virtual bool initialize() override
    {
      sample0Item::initializeClass(this);
      sample1Item::initializeClass(this);
      return true;
    }
  };
}

CNOID_IMPLEMENT_PLUGIN_ENTRY(ik_constraint2_body_contact_sample::IKConstraintBodyContactSamplePlugin)
