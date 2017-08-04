#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "hstcalversion.h"

char * getVersionInfo(char ** buffer)
{
    if (!buffer)
        return NULL;
    if (*buffer)
        assert(0); // Incorrect usage - NULL ptr must be passed in.

    const char * format = "Running %s-%s, git branch: %s HEAD @: %s";
    size_t length = strlen(format) + strlen(APPNAME) + strlen(VERSION) + strlen(BRANCH) + strlen(COMMIT); // NOTE: this doesn't require an explicit +1 for '\0' as it is 8 larger than needed from '%s'.
    *buffer = malloc(length*sizeof(char));
    if (!*buffer)
        return NULL;
    sprintf(*buffer, format, APPNAME, VERSION, BRANCH, COMMIT);
    return *buffer;
}

//Move this to common hstcal/lib/trlbuf.c when finished #191
void trlVersion(void)
{
    char * versionText = NULL;
    getVersionInfo(&versionText);
    printf("%s\n", versionText);//trlmessage(versionText); //waiting on #191
    if (versionText)
    {
        free(versionText);
        versionText = NULL;
    }
}
