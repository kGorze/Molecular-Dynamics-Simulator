#include <stdio.h>
#include <string.h>

void GetNameList(const char *fd) {
    FILE *file = fopen(fd, "r");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    
    char line[256];
    char *key, *value;
    char *pattern = "initUcell";
    while (fgets(line, sizeof(line), file) != NULL) {
        if (strncmp(line, pattern, strlen(pattern)) == 0) {
            key = strtok(line + strlen(pattern), " \t\n");
            value = strtok(NULL, " \t\n");
            // Matrix of molecular unit cells
            // _mdsim_globals[k] = _namelist_converter[k](nx, ny);
            printf("%s %s\n", key, value); // Placeholder for processing the key-value pair
        } else {
            key = strtok(line, " \t\n");
            value = strtok(NULL, " \t\n");
            // _mdsim_globals[k] = _namelist_converter[k](v);
            printf("%s %s\n", key, value); // Placeholder for processing the key-value pair
        }
    }
    
    fclose(file);
}

int main() {
    GetNameList("data.in"); // Replace "filename.txt" with the actual filename
    return 0;
}
