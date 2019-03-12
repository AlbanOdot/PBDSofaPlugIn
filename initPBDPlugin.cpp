#include "initPBDPlugin.hpp"

extern "C" {
    void initExternalModule()
    {
        static bool first = true;
        if (first)
        {
            first = false;
        }
    }

    const char* getModuleName()
    {
        return "PBD AlphaV";
    }

    const char* getModuleVersion()
    {
        return "0.0.1";
    }

    const char* getModuleLicense()
    {
        return "LGPL";
    }

    const char* getModuleDescription()
    {
        return "Basic PBD implementation";
    }

    const char* getModuleComponentList()
    {
        // Comma-separated list of the components in this plugin, empty for now
        return "PBDAnimationLoop";
    }
}
