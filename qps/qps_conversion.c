#include "qpalm.h"
#include "constants.h"
#include "global_opts.h"
#include "cholmod.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define BUFFER_LEN 9 

void c_strcpy(char dest[], const char source[], size_t len) {
    size_t i;
    for(i = 0; i < len && (dest[i] = source[i]) != '\0'; i++);
    dest[len-1] = '\0';
}

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

    char line[100], command[20], name[50], NLGE[1], buf1[9], buf1_copy[9], buf2[9], buf2_copy[9], buf3[9], buf3_copy[9], temp1[BUFFER_LEN+1], temp2[BUFFER_LEN+1];
    size_t len = 0;

    fgets(line, 100, fp_in);
    fputs(line, fp_out);
    char next_char;
    next_char = fgetc(fp_in);
    // get_next_command(command, next_char, fp_in);
    // command[0] = next_char;

    // fgets(line, 100, fp_in);
    // sscanf(line, "%s", &command[1]);
    // fputs(command, fp_out);
    // fputs("\n", fp_out);

    // next_char = fgetc(fp_in);

    while (get_next_command(command, next_char, fp_in)) {
      fputs(command, fp_out);
      fputs("\n", fp_out);
      next_char = fgetc(fp_in);
      // printf("Command: %s\n", command);
      if (!strcmp(command, "ROWS")) {
        while (next_char == ' ') {
            fgets(line, 100, fp_in);
            c_strcpy(buf1, &line[3], BUFFER_LEN);
            remove_spaces(buf1_copy, buf1);
            line[3] = '\0';
            strcat(line, buf1_copy);
            fputs(" ", fp_out);
            fputs(line, fp_out);
            fputs("\n", fp_out);
            // printf("Line:%s\n", line);
            next_char = fgetc(fp_in);
        }
      } else if (!strcmp(command, "COLUMNS")) {
          while (next_char == ' ') {
            fgets(line, 100, fp_in);
            len = strlen(line);
            c_strcpy(buf1, &line[3], BUFFER_LEN);
            remove_spaces(buf1_copy, buf1);
            c_strcpy(buf2, &line[13], BUFFER_LEN);
            remove_spaces(buf2_copy, buf2);
            c_strcpy(temp1, &line[27], BUFFER_LEN+1);
            if (len > 40) {
              c_strcpy(buf3, &line[38], BUFFER_LEN);
              remove_spaces(buf3_copy, buf3);
              c_strcpy(temp2, &line[51], BUFFER_LEN+1);
            }

            line[0] = '\0';
            strcat(line, buf1_copy);
            strcat(line, "  ");
            strcat(line, buf2_copy);
            strcat(line, "    ");
            strcat(line, temp1);
            if (len > 40) {
              strcat(line, "    ");
              strcat(line, buf3_copy);
              strcat(line, "    ");
              strcat(line, temp2);
            }
            fputs("    ", fp_out);
            fputs(line, fp_out);
            fputs("\n", fp_out);
            next_char = fgetc(fp_in);
        }
      } else if (!strcmp(command, "RHS")) {
        while (next_char == ' ') {
            fgets(line, 100, fp_in);
            len = strlen(line);
            // printf("line: %s\n", line);
            c_strcpy(buf1, &line[3], BUFFER_LEN);
            remove_spaces(buf1_copy, buf1);
            c_strcpy(buf2, &line[13], BUFFER_LEN);
            remove_spaces(buf2_copy, buf2);
            c_strcpy(temp1, &line[27], BUFFER_LEN+1);
            if (len > 40) {
              c_strcpy(buf3, &line[38], BUFFER_LEN);
              remove_spaces(buf3_copy, buf3);
              // printf("buf3: %s, buf3_copy: %s\n", buf3, buf3_copy);
              c_strcpy(temp2, &line[51], BUFFER_LEN+1);
            }

            line[0] = '\0';
            strcat(line, buf1_copy);
            strcat(line, "  ");
            strcat(line, buf2_copy);
            strcat(line, "    ");
            strcat(line, temp1);
            if (len > 40) {
              strcat(line, "    ");
              strcat(line, buf3_copy);
              strcat(line, "    ");
              strcat(line, temp2);
            }
            fputs("    ", fp_out);
            fputs(line, fp_out);
            fputs("\n", fp_out);
            next_char = fgetc(fp_in);
        }
      } else if (!strcmp(command, "RANGES")) {
        while (next_char == ' ') {
            fgets(line, 100, fp_in);
            len = strlen(line);
            // printf("line: %s\n", line);
            c_strcpy(buf1, &line[3], BUFFER_LEN);
            remove_spaces(buf1_copy, buf1);
            c_strcpy(buf2, &line[13], BUFFER_LEN);
            remove_spaces(buf2_copy, buf2);
            c_strcpy(temp1, &line[27], BUFFER_LEN+1);
            if (len > 40) {
              c_strcpy(buf3, &line[38], BUFFER_LEN);
              remove_spaces(buf3_copy, buf3);
              // printf("buf3: %s, buf3_copy: %s\n", buf3, buf3_copy);
              c_strcpy(temp2, &line[51], BUFFER_LEN+1);
            }

            line[0] = '\0';
            strcat(line, buf1_copy);
            strcat(line, "  ");
            strcat(line, buf2_copy);
            strcat(line, "    ");
            strcat(line, temp1);
            if (len > 40) {
              // printf("Len: %u, line: %s\n", len, line);
              strcat(line, "    ");
              strcat(line, buf3_copy);
              strcat(line, "    ");
              strcat(line, temp2);
            }
            fputs("    ", fp_out);
            fputs(line, fp_out);
            fputs("\n", fp_out);
            next_char = fgetc(fp_in);
        }
      } else if (!strcmp(command, "BOUNDS")) {
        while (next_char == ' ') {
            fgets(line, 100, fp_in);
            len = strlen(line);
            // printf("line: %s\n", line);
            c_strcpy(buf1, &line[3], BUFFER_LEN);
            remove_spaces(buf1_copy, buf1);
            c_strcpy(buf2, &line[13], BUFFER_LEN);
            remove_spaces(buf2_copy, buf2);
            c_strcpy(temp1, &line[27], BUFFER_LEN+1);
            // if (len > 40) {
            //   c_strcpy(buf3, &line[38], BUFFER_LEN);
            //   remove_spaces(buf3_copy, buf3);
            //   // printf("buf3: %s, buf3_copy: %s\n", buf3, buf3_copy);
            //   c_strcpy(temp2, &line[51], BUFFER_LEN+1);
            // }

            line[3] = '\0';
            strcat(line, buf1_copy);
            strcat(line, "  ");
            strcat(line, buf2_copy);
            strcat(line, "    ");
            strcat(line, temp1);
            // if (len > 40) {
            //   // printf("Len: %u, line: %s\n", len, line);
            //   strcat(line, "    ");
            //   strcat(line, buf3_copy);
            //   strcat(line, "    ");
            //   strcat(line, temp2);
            // }
            // fputs("    ", fp_out);
            fputs(" ", fp_out);
            fputs(line, fp_out);
            fputs("\n", fp_out);
            next_char = fgetc(fp_in);
        }
      } else {
        break;
      }
    }
    fputs("ENDATA\n", fp_out);



    // while (next_char == ' ') {
    //     fgets(line, 100, fp_in);
    //     c_strcpy(buf1, &line[3], BUFFER_LEN);
    //     remove_spaces(buf1_copy, buf1);
    //     line[3] = '\0';
    //     strcat(line, buf1_copy);
    //     fputs(" ", fp_out);
    //     fputs(line, fp_out);
    //     fputs("\n", fp_out);
    //     next_char = fgetc(fp_in);
    // }

    // get_next_command(command, next_char, fp_in);
    // fputs(command, fp_out);
    // fputs("\n", fp_out);
    // next_char = fgetc(fp_in);
    // printf("Reading %s\n", command);
    // while (next_char == ' ') {
    //     fgets(line, 100, fp_in);
    //     len = strlen(line);
    //     c_strcpy(buf1, &line[3], BUFFER_LEN);
    //     remove_spaces(buf1_copy, buf1);
    //     c_strcpy(buf2, &line[13], BUFFER_LEN);
    //     remove_spaces(buf2_copy, buf2);
    //     c_strcpy(temp1, &line[29], BUFFER_LEN);
    //     if (len > 40) {
    //       c_strcpy(buf3, &line[38], BUFFER_LEN);
    //       remove_spaces(buf3_copy, buf3);
    //       c_strcpy(temp2, &line[54], BUFFER_LEN);
    //     }

    //     line[0] = '\0';
    //     strcat(line, buf1_copy);
    //     strcat(line, "  ");
    //     strcat(line, buf2_copy);
    //     strcat(line, "    ");
    //     strcat(line, temp1);
    //     if (len > 40) {
    //       strcat(line, buf3_copy);
    //       strcat(line, "    ");
    //       strcat(line, temp2);
    //     }
    //     fputs("    ", fp_out);
    //     fputs(line, fp_out);

    //     fputs("\n", fp_out);
    //     next_char = fgetc(fp_in);
    // }
    // get_next_command(command, next_char, fp_in);
    // fputs(command, fp_out);
    // fputs("\n", fp_out);
    // next_char = fgetc(fp_in);
    // printf("Reading %s\n", command);

    // while (next_char == ' ') {


    //   next_char = fgetc(fp_in);
    // }


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

    c_free(new_filename);
    return 0;

}