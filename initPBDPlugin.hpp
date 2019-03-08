#ifndef INITMYPLUGIN_H
#define INITMYPLUGIN_H

#include <sofa/helper/system/config.h>

#ifdef SOFA_BUILD_MYPLUGIN
#define SOFA_MyPlugin_API SOFA_EXPORT_DYNAMIC_LIBRARY
#else
#define SOFA_MyPlugin_API SOFA_IMPORT_DYNAMIC_LIBRARY
#endif

/** mainpage
	NO DOC HEHEHEHEHE *VOICE BECOMING DARKER* MOUHAHAHAHAHAHHAA *kuff kuff*
	Sorry dear sir no doc
 */

#endif
