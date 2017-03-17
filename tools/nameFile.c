// Utility for converting parts of a file name into a full file name
// Kelly McQuighan
// 8/4/13
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nameFile.h"

void nameFile ( char* folder, char* prefix, char* name, char* suffix, char* filename )
{
int folder_len, prefix_len, name_len, suffix_len, i;

    folder_len = strlen(folder);
    prefix_len = strlen(prefix);
    name_len   = strlen(name);
    suffix_len = strlen(suffix);

    for (i=0; i<folder_len; i++) {
      filename[i] = folder[i];
    }
    filename[folder_len] = '/';
    for (i=0; i<prefix_len; i++) {
      filename[i+folder_len+1]= prefix[i];
    }
    for (i=0; i<name_len; i++) {
      filename[i+folder_len+prefix_len+1]= name[i];
    }
    filename[folder_len+prefix_len+name_len+1] = '.';
    for (i=0; i<suffix_len; i++) {
      filename[i+folder_len+prefix_len+name_len+2]= suffix[i];
    }
    filename[folder_len+prefix_len+name_len+suffix_len+2] = '\0';
}
