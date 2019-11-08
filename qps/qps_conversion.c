#include "qpalm.h"
#include "constants.h"
#include "global_opts.h"
#include "cholmod.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
void remove_spaces (char* restrict str_trimmed, const char* restrict str_untrimmed)
{
  while (*str_untrimmed != '\0')
  {
    if(!isspace(*str_untrimmed))
    {
      *str_trimmed = *str_untrimmed;
      str_trimmed++;
    }
    str_untrimmed++;
  }
  *str_trimmed = '\0';
}

int get_next_command(char* command, char next_char, FILE* fp) {
    command[0] = next_char;
    char line[100];
    fgets(line, 100, fp);
    sscanf(line, "%s", &command[1]);

    return strcmp(command, "ENDATA");
}

void convert_to_new_format(const char* filename, char* new_filename) {
    FILE* fp_in, *fp_out;

    fp_in = fopen(filename, "r");
    if(fp_in == NULL) {
        fprintf(stderr, "Could not open file %s\n", filename);
        return 1;
    }
    fp_out = fopen(new_filename, "w");
    if(fp_out == NULL) {
        fprintf(stderr, "Could not open file %s\n", filename);
        return 1;
    }

    char line[100], command[20], name[50], NLGE[1], buf[20], buf_copy[20];
    fgets(line, 100, fp_in);
    fputs(line, fp_out);
    char next_char;
    next_char = fgetc(fp_in);
    command[0] = next_char;

    fgets(line, 100, fp_in);
    sscanf(line, "%s", &command[1]);
    fputs(command, fp_out);
    fputs("\n", fp_out);

    next_char = fgetc(fp_in);
    while (next_char == ' ') {
        fgets(line, 100, fp_in);
        strcpy(buf, &line[5]);
        // sscanf(line, "%s %s", NLGE, buf);
        remove_spaces(buf_copy, buf);
        line[5] = '\0';
        strcat(line, buf_copy);
        fputs(" ", fp_out);
        fputs(line, fp_out);
        // fputs(" ", fp_out);
        // fputs(NLGE, fp_out);
        // fputs("  ", fp_out);
        // fputs(buf_copy, fp_out);
        fputs("\n", fp_out);
        next_char = fgetc(fp_in);
    }


}

int main(int argc, char *argv[]) {
    char* filename = argv[1];
    char* new_filename = c_malloc(strlen(filename)+5+1);
    new_filename[0] = '\0';
    strcat(new_filename, filename);
    new_filename[strlen(filename)-4] = '\0';
    strcat(new_filename, "_copy");
    strcat(new_filename, ".qps");
    convert_to_new_format(filename, new_filename);  
    return 0;

}