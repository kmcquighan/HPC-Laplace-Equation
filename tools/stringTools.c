// Conversions to character arrays because these and sstring don't work
// on my computer.
// Kelly McQuighan
// last modified 8/1/13

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include "stringTools.h"

//*************Int to Char*********************
void itoa(int num, char* buffer) {
int l, i, neg=0, base=10;
  if ( num < 0 ) { num *= -1; neg = 1; }
  if ( num == 0 ) {
    l = 1;
    buffer[0] = '0'; buffer[1] = '\0';
    return;
  } 

  l = log10(num) + 1; // Length of number like if num=123456 then l=6.
  if ( neg == 1 ) { l++; buffer[0] = '-'; }
  i = l - 1;

  while(num>0 && i>=0)
  {
    buffer[i--]=(char)(num%base+48);
    num/=base;
  }
  buffer[l]='\0';
//  fprintf(stdout, "The number is %s.\n",buffer);

return;
}

//***********Int to char, returns int length*********
int itoa_len(int num, char* buffer) {
int l, i, neg = 0, base=10;
  if ( num < 0 ) { num *= -1; neg = 1; }
  if ( num == 0 ) {
    l = 1;
    buffer[0] = '0'; buffer[1] = '\0';
    return l;
  } 
  l = log10(num) + 1; // Length of number like if num=123456 then l=6.
  if ( neg == 1 ) { l++; buffer[0] = '-'; }
  i = l - 1;

  while(num>0 && i>=0)
  {
    buffer[i--]=(char)(num%base+48);
    num/=base;
  }
  buffer[l]='\0';
//  fprintf(stdout, "The number is %s.\n",buffer);

return l;
}

//*******Double to char***************************
void  dtoa(double x, char* buffer, char* format, int numDecimal) {
int neg = 0, xtemp, ltemp, l, power_len, headlen, st=1, i, power;
int diff, base=10;
char temp[33], powerbuff[33];

  if (numDecimal < 0 ) {fprintf(stdout, "ERROR: number of decimal places must be at least zero.\n"); return;}
  if ( x < 0 ) { x *= -1; neg = 1; }
  if ( x == 0 ) {
    l = 3+numDecimal;
    buffer[0] = ' '; buffer[1] = '0'; buffer[2] = '.';
    for (i=3; i<l; i++) buffer[i] = '0';
    if ( (format[0]=='e') || (format[0]=='E') ) {
      buffer[l] = 'e'; buffer[l+1]='+'; buffer[l+2]='0'; buffer[l+3]='1'; buffer[l+4]='\0';
    } else if ( (format[0]=='f') || (format[0]=='F') ) {
      buffer[l] = '\0';
    } else {
      fprintf(stdout, "ERROR: format choices are 'e' for scientific notation or 'f' for decimal notation\n");
    }
    return;
  }

  if ( (format[0]=='e') || (format[0]=='E') ) {
    if ( floor(x) > 0 ) {
      xtemp = floor(x);
      power = itoa_len(xtemp, temp) - 1;
      diff = (numDecimal - power);
      if (diff > 0) {
        xtemp = floor(x*pow(base, diff));
        itoa(xtemp, temp);
      }
      if ( neg == 1 ) buffer[0] = '-';
      else            buffer[0] = ' ';
      buffer[1] = temp[0]; buffer[2] = '.';
      for (i=1; i<=numDecimal; i++) buffer[i+2] = temp[i];
      buffer[numDecimal+3] = 'e';
      buffer[numDecimal+4] = '+';
      power_len = itoa_len(power, powerbuff);
      if (power_len==1) { 
        buffer[numDecimal+5] = '0'; 
        buffer[numDecimal+6] = powerbuff[0];
        buffer[numDecimal+7] = '\0';
      } else {
        for (i=0; i<power_len; i++) buffer[numDecimal+5+i] = powerbuff[i];
        buffer[numDecimal+5+power_len] = '\0';
      }
    } else {
      xtemp = floor(x);
      power=0;
      while (xtemp == 0) {
        power++;
        xtemp = floor(x*pow(base,power));
      }
      xtemp = floor(x*pow(base,power+numDecimal+1));
      itoa_len(xtemp, temp);
      if ( neg == 1 ) buffer[0] = '-';
      else            buffer[0] = ' ';
      buffer[1] = temp[0]; buffer[2] = '.';
      for (i=1; i<=numDecimal; i++) buffer[i+2] = temp[i];
      buffer[numDecimal+3] = 'e';
      buffer[numDecimal+4] = '-';
      power_len = itoa_len(power, powerbuff);
      if (power_len==1) { 
        buffer[numDecimal+5] = '0'; 
        buffer[numDecimal+6] = powerbuff[0];
        buffer[numDecimal+7] = '\0';
      } else {
        for (i=0; i<power_len; i++) buffer[numDecimal+5+i] = powerbuff[i];
        buffer[numDecimal+5+power_len] = '\0';
      }
    }

  } else if ( (format[0]=='f') || (format[0]=='F') ) {
    xtemp = floor(x*pow(base, numDecimal));
    ltemp = itoa_len(xtemp, temp);
    headlen = ltemp - numDecimal;
    l = ltemp + 2;
    if ( neg == 1 ) {
      buffer[0] = '-';
      if ( headlen <= 0 ) { l++; st++; buffer[1] = '0'; buffer[2]='.' ;}
    } else {
      buffer[0] = ' ';
      if ( headlen <= 0 ) { l++; st++; buffer[1] = '0'; buffer[2]='.'; }
    }
 //   fprintf(stdout, "%d\n", l);
    for ( i=0; i<headlen; i++ ) buffer[i+st] = temp[i];
    if (headlen > 0 ) buffer[headlen+st] = '.';
    while ( headlen < 0 ) {
      headlen++; st++; l++;
      buffer[st] = '0';
    }
    for ( i=headlen; i<ltemp; i++ ) buffer[i+st+1] = temp[i];
    buffer[l]='\0';

  } else {
    fprintf(stdout, "ERROR: format choices are 'e' for scientific notation or 'f' for decimal notation\n");
  }
return;
}

//*********Double to char returns length*****************
int  dtoa_len(double x, char* buffer, char* format, int numDecimal) {
int neg = 0, xtemp, ltemp, l, power_len, headlen, st=1, i, power;
int diff, base=10;
char temp[33], powerbuff[33];

  if (numDecimal < 0 ) {fprintf(stdout, "ERROR: number of decimal places must be at least zero.\n"); return 0;}
  if ( x < 0 ) { x *= -1; neg = 1; }
  if ( x == 0 ) {
    l = 3+numDecimal;
    buffer[0] = ' '; buffer[1] = '0'; buffer[2] = '.';
    for (i=3; i<l; i++) buffer[i] = '0';
    if ( (format[0]=='e') || (format[0]=='E') ) {
      buffer[l] = 'e'; buffer[l+1]='+'; buffer[l+2]='0'; buffer[l+3]='1'; buffer[l+4]='\0';
      return l+4;
    } else if ( (format[0]=='f') || (format[0]=='F') ) {
      buffer[l] = '\0';
      return l;
    } else {
      fprintf(stdout, "ERROR: format choices are 'e' for scientific notation or 'f' for decimal notation\n");
      return 0;
    }
  } 

  if ( (format[0]=='e') || (format[0]=='E') ) {
    if ( floor(x) > 0 ) {
      xtemp = floor(x);
      power = itoa_len(xtemp, temp) - 1;
      diff = (numDecimal - power);
      if (diff > 0) {
        xtemp = floor(x*pow(base, diff));
        itoa(xtemp, temp);
      }
      if ( neg == 1 ) buffer[0] = '-';
      else            buffer[0] = ' ';
      buffer[1] = temp[0]; buffer[2] = '.';
      for (i=1; i<=numDecimal; i++) buffer[i+2] = temp[i];
      buffer[numDecimal+3] = 'e';
      buffer[numDecimal+4] = '+';
      power_len = itoa_len(power, powerbuff);
      if (power_len==1) { 
        buffer[numDecimal+5] = '0'; 
        buffer[numDecimal+6] = powerbuff[0];
        buffer[numDecimal+7] = '\0';
        l = numDecimal+7;
      } else {
        for (i=0; i<power_len; i++) buffer[numDecimal+5+i] = powerbuff[i];
        buffer[numDecimal+5+power_len] = '\0';
        l=numDecimal+5+power_len;
      }
    } else {
      xtemp = floor(x);
      power=0;
      while (xtemp == 0) {
        power++;
        xtemp = floor(x*pow(base,power));
      }
      xtemp = floor(x*pow(base,power+numDecimal+1));
      itoa_len(xtemp, temp);
      if ( neg == 1 ) buffer[0] = '-';
      else            buffer[0] = ' ';
      buffer[1] = temp[0]; buffer[2] = '.';
      for (i=1; i<=numDecimal; i++) buffer[i+2] = temp[i];
      buffer[numDecimal+3] = 'e';
      buffer[numDecimal+4] = '-';
      power_len = itoa_len(power, powerbuff);
      if (power_len==1) { 
        buffer[numDecimal+5] = '0'; 
        buffer[numDecimal+6] = powerbuff[0];
        buffer[numDecimal+7] = '\0';
        l = numDecimal+7;
      } else {
        for (i=0; i<power_len; i++) buffer[numDecimal+5+i] = powerbuff[i];
        buffer[numDecimal+5+power_len] = '\0';
        l = numDecimal+5+power_len;
      }
    }

  } else if ( (format[0]=='f') || (format[0]=='F') ) {
    xtemp = floor(x*pow(base, numDecimal));
    ltemp = itoa_len(xtemp, temp);
    headlen = ltemp - numDecimal;
    l = ltemp + 2;
    if ( neg == 1 ) {
      buffer[0] = '-';
      if ( headlen <= 0 ) { l++; st++; buffer[1] = '0'; buffer[2]='.' ;}
    } else {
      buffer[0] = ' ';
      if ( headlen <= 0 ) { l++; st++; buffer[1] = '0'; buffer[2]='.'; }
    }
 //   fprintf(stdout, "%d\n", l);
    for ( i=0; i<headlen; i++ ) buffer[i+st] = temp[i];
    if (headlen > 0 ) buffer[headlen+st] = '.';
    while ( headlen < 0 ) {
      headlen++; st++; l++;
      buffer[st] = '0';
    }
    for ( i=headlen; i<ltemp; i++ ) buffer[i+st+1] = temp[i];
    buffer[l]='\0';

  } else {
    fprintf(stdout, "ERROR: format choices are 'e' for scientific notation or 'f' for decimal notation\n");
    l=0;
  }

return l;
}
//*******Single to char***************************
void  stoa(float x, char* buffer, char* format, int numDecimal) {
int neg = 0, xtemp, ltemp, l, power_len, headlen, st=1, i, power;
int diff, base=10;
char temp[33], powerbuff[33];

  if (numDecimal < 0 ) {fprintf(stdout, "ERROR: number of decimal places must be at least zero.\n"); return;}
  if ( x < 0 ) { x *= -1; neg = 1; }
  if ( x == 0 ) {
    l = 3+numDecimal;
    buffer[0] = ' '; buffer[1] = '0'; buffer[2] = '.';
    for (i=3; i<l; i++) buffer[i] = '0';
    if ( (format[0]=='e') || (format[0]=='E') ) {
      buffer[l] = 'e'; buffer[l+1]='+'; buffer[l+2]='0'; buffer[l+3]='1'; buffer[l+4]='\0';
    } else if ( (format[0]=='f') || (format[0]=='F') ) {
      buffer[l] = '\0';
    } else {
      fprintf(stdout, "ERROR: format choices are 'e' for scientific notation or 'f' for decimal notation\n");
    }
    return;
  }

  if ( (format[0]=='e') || (format[0]=='E') ) {
    if ( floor(x) > 0 ) {
      xtemp = floor(x);
      power = itoa_len(xtemp, temp) - 1;
      diff = (numDecimal - power);
      if (diff > 0) {
        xtemp = floor(x*pow(base, diff));
        itoa(xtemp, temp);
      }
      if ( neg == 1 ) buffer[0] = '-';
      else            buffer[0] = ' ';
      buffer[1] = temp[0]; buffer[2] = '.';
      for (i=1; i<=numDecimal; i++) buffer[i+2] = temp[i];
      buffer[numDecimal+3] = 'e';
      buffer[numDecimal+4] = '+';
      power_len = itoa_len(power, powerbuff);
      if (power_len==1) { 
        buffer[numDecimal+5] = '0'; 
        buffer[numDecimal+6] = powerbuff[0];
        buffer[numDecimal+7] = '\0';
      } else {
        for (i=0; i<power_len; i++) buffer[numDecimal+5+i] = powerbuff[i];
        buffer[numDecimal+5+power_len] = '\0';
      }
    } else {
      xtemp = floor(x);
      power=0;
      while (xtemp == 0) {
        power++;
        xtemp = floor(x*pow(base,power));
      }
      xtemp = floor(x*pow(base,power+numDecimal+1));
      itoa_len(xtemp, temp);
      if ( neg == 1 ) buffer[0] = '-';
      else            buffer[0] = ' ';
      buffer[1] = temp[0]; buffer[2] = '.';
      for (i=1; i<=numDecimal; i++) buffer[i+2] = temp[i];
      buffer[numDecimal+3] = 'e';
      buffer[numDecimal+4] = '-';
      power_len = itoa_len(power, powerbuff);
      if (power_len==1) { 
        buffer[numDecimal+5] = '0'; 
        buffer[numDecimal+6] = powerbuff[0];
        buffer[numDecimal+7] = '\0';
      } else {
        for (i=0; i<power_len; i++) buffer[numDecimal+5+i] = powerbuff[i];
        buffer[numDecimal+5+power_len] = '\0';
      }
    }

  } else if ( (format[0]=='f') || (format[0]=='F') ) {
    xtemp = floor(x*pow(base, numDecimal));
    ltemp = itoa_len(xtemp, temp);
    headlen = ltemp - numDecimal;
    l = ltemp + 2;
    if ( neg == 1 ) {
      buffer[0] = '-';
      if ( headlen <= 0 ) { l++; st++; buffer[1] = '0'; buffer[2]='.' ;}
    } else {
      buffer[0] = ' ';
      if ( headlen <= 0 ) { l++; st++; buffer[1] = '0'; buffer[2]='.'; }
    }
 //   fprintf(stdout, "%d\n", l);
    for ( i=0; i<headlen; i++ ) buffer[i+st] = temp[i];
    if (headlen > 0 ) buffer[headlen+st] = '.';
    while ( headlen < 0 ) {
      headlen++; st++; l++;
      buffer[st] = '0';
    }
    for ( i=headlen; i<ltemp; i++ ) buffer[i+st+1] = temp[i];
    buffer[l]='\0';

  } else {
    fprintf(stdout, "ERROR: format choices are 'e' for scientific notation or 'f' for decimal notation\n");
  }
return;
}

//*********Single to char returns length*****************
int  stoa_len(float x, char* buffer, char* format, int numDecimal) {
int neg = 0, xtemp, ltemp, l, power_len, headlen, st=1, i, power;
int diff, base=10;
char temp[33], powerbuff[33];

  if (numDecimal < 0 ) {fprintf(stdout, "ERROR: number of decimal places must be at least zero.\n"); return 0;}
  if ( x < 0 ) { x *= -1; neg = 1; }
  if ( x == 0 ) {
    l = 3+numDecimal;
    buffer[0] = ' '; buffer[1] = '0'; buffer[2] = '.';
    for (i=3; i<l; i++) buffer[i] = '0';
    if ( (format[0]=='e') || (format[0]=='E') ) {
      buffer[l] = 'e'; buffer[l+1]='+'; buffer[l+2]='0'; buffer[l+3]='1'; buffer[l+4]='\0';
      return l+4;
    } else if ( (format[0]=='f') || (format[0]=='F') ) {
      buffer[l] = '\0';
      return l;
    } else {
      fprintf(stdout, "ERROR: format choices are 'e' for scientific notation or 'f' for decimal notation\n");
      return 0;
    }
  } 

  if ( (format[0]=='e') || (format[0]=='E') ) {
    if ( floor(x) > 0 ) {
      xtemp = floor(x);
      power = itoa_len(xtemp, temp) - 1;
      diff = (numDecimal - power);
      if (diff > 0) {
        xtemp = floor(x*pow(base, diff));
        itoa(xtemp, temp);
      }
      if ( neg == 1 ) buffer[0] = '-';
      else            buffer[0] = ' ';
      buffer[1] = temp[0]; buffer[2] = '.';
      for (i=1; i<=numDecimal; i++) buffer[i+2] = temp[i];
      buffer[numDecimal+3] = 'e';
      buffer[numDecimal+4] = '+';
      power_len = itoa_len(power, powerbuff);
      if (power_len==1) { 
        buffer[numDecimal+5] = '0'; 
        buffer[numDecimal+6] = powerbuff[0];
        buffer[numDecimal+7] = '\0';
        l = numDecimal+7;
      } else {
        for (i=0; i<power_len; i++) buffer[numDecimal+5+i] = powerbuff[i];
        buffer[numDecimal+5+power_len] = '\0';
        l=numDecimal+5+power_len;
      }
    } else {
      xtemp = floor(x);
      power=0;
      while (xtemp == 0) {
        power++;
        xtemp = floor(x*pow(base,power));
      }
      xtemp = floor(x*pow(base,power+numDecimal+1));
      itoa_len(xtemp, temp);
      if ( neg == 1 ) buffer[0] = '-';
      else            buffer[0] = ' ';
      buffer[1] = temp[0]; buffer[2] = '.';
      for (i=1; i<=numDecimal; i++) buffer[i+2] = temp[i];
      buffer[numDecimal+3] = 'e';
      buffer[numDecimal+4] = '-';
      power_len = itoa_len(power, powerbuff);
      if (power_len==1) { 
        buffer[numDecimal+5] = '0'; 
        buffer[numDecimal+6] = powerbuff[0];
        buffer[numDecimal+7] = '\0';
        l = numDecimal+7;
      } else {
        for (i=0; i<power_len; i++) buffer[numDecimal+5+i] = powerbuff[i];
        buffer[numDecimal+5+power_len] = '\0';
        l = numDecimal+5+power_len;
      }
    }

  } else if ( (format[0]=='f') || (format[0]=='F') ) {
    xtemp = floor(x*pow(base, numDecimal));
    ltemp = itoa_len(xtemp, temp);
    headlen = ltemp - numDecimal;
    l = ltemp + 2;
    if ( neg == 1 ) {
      buffer[0] = '-';
      if ( headlen <= 0 ) { l++; st++; buffer[1] = '0'; buffer[2]='.' ;}
    } else {
      buffer[0] = ' ';
      if ( headlen <= 0 ) { l++; st++; buffer[1] = '0'; buffer[2]='.'; }
    }
 //   fprintf(stdout, "%d\n", l);
    for ( i=0; i<headlen; i++ ) buffer[i+st] = temp[i];
    if (headlen > 0 ) buffer[headlen+st] = '.';
    while ( headlen < 0 ) {
      headlen++; st++; l++;
      buffer[st] = '0';
    }
    for ( i=headlen; i<ltemp; i++ ) buffer[i+st+1] = temp[i];
    buffer[l]='\0';

  } else {
    fprintf(stdout, "ERROR: format choices are 'e' for scientific notation or 'f' for decimal notation\n");
    l=0;
  }

return l;
}
