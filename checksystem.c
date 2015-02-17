# include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define memfile "/proc/meminfo"

typedef struct datamem{
    unsigned long int MemTotal, MemFree, Cached, Buffers, MemAvail, MemUsed;
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
    in->MemAvail=in->MemFree+in->Cached+in->Buffers;
    in->MemUsed=in->MemTotal-in->MemAvail;
    fclose(doc);
    return 0;
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
        printf("MemAvail %lu\n",data.MemAvail);
        printf("MemUsed %lu\n",data.MemUsed);
        }
	return 0;
}

