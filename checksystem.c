/******************************************************
 * This file is part of the Codespeckle++ distribution Copyright (c) 2015 Jimmy
 * Aguilar Mena.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define memfile "/proc/meminfo"

typedef struct datamem{
    unsigned long int MemTotal, MemFree, Cached, Buffers;
    unsigned long int MemAvail(){return MemFree+Cached+Buffers;}
    unsigned long int MemUsed(){return MemTotal-MemFree-Cached-Buffers;}
    } datamem;

int getmem(datamem* in){
	char index[20], unit[2], temp[100];
    unsigned long int value;
    FILE *doc = fopen(memfile,"r");
    if (doc==NULL) return -1;
    while (!feof(doc)){
        fgets(temp,100,doc);
        sscanf(temp,"%s %lu %s",index,&value,unit);
        if (strncmp(index,"MemTotal:",9)==0) in->MemTotal=value;
        else if (strncmp(index,"MemFree:",8)==0) in->MemFree=value;
        else if (strncmp(index,"Cached:",7)==0) in->Cached=value;
        else if (strncmp(index,"Buffers:",8)==0) in->Buffers=value;
        }
    fclose(doc);
    return 0;
    }

int getcpu(){
    
    }

int main(int argc, char **argv){
    datamem data;
    if (getmem(&data)==-1){
        printf("Error opening memory data file: %s", memfile);
        exit(-1);
        }
    else{
        printf("MemTotal %lu\n",data.MemTotal);
        printf("MemFree %lu\n",data.MemFree);
        printf("Cached %lu\n",data.Cached);
        printf("Buffers %lu\n",data.Buffers);
        printf("MemAvail %lu\n",data.MemAvail());
        printf("MemUsed %lu\n",data.MemUsed());
        }
	return 0;
}

