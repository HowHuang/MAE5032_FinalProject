#include <stdio.h>
#include <stdlib.h>
int main()
{   
    char buf[1000];
    int i = 99;
    sprintf(buf,"%08d",i);
    printf("%s\n",buf);
    printf("atoi(00000099)=%d\n",atoi(buf));
}