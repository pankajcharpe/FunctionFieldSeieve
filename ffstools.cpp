/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include"ffstools.h"



//Set pol to t^i
void GF2X_set_ti(GF2X& pol,long i)
{
   GF2X f;
   set(f);
   f=LeftShift(f,i);
   pol=f;
   f.kill();
}


void GF2X_submul(GF2X& f,GF2X& g,GF2X& h)
{
  GF2X t;
  t=g*h;
  f=f-t;
  t.kill();
}
 
//GF2X to hex
void GF2X_to_hex(GF2X f)
{

  ofstream outfile;
   outfile.open("date.txt",ios::app);
   
  char hex[100];
  long i=1, j, temp;
  long d=deg(f);
  long index=0;
  f=reverse(f);
  ZZ decimal=ZZ(0);
  while(index<=d)
  {
    ZZ tmp=(power2_ZZ(index)*(ZZ)rep(f[index]));
    decimal=decimal+tmp;
    index++;
  }
   while(decimal!=ZZ(0))
	{
		temp = decimal % long(16);
	
		if( temp < 10)
		{
			temp = temp + 48;
		}
		else
		{
			temp = temp + 55;
		}
		hex[i++]= temp;
		decimal = decimal / long(16);
	}
	
	for(j=i-1 ;j>0;j--)
	{
		outfile<<hex[j];
	}
	 
  
} 
 
uint64_t power(long i)
{
  uint64_t num=1;
  for(long j=0;j<i;j++)
   num=num*2;
  return num;
}

//Conversion from GF2X to uint64_t

uint64_t GF2X_to_uint64_t(GF2X& f)
{
  uint64_t num=0;
  long d=deg(f);
  long i=0;
  uint64_t base=0;
  GF2 ith_bit;
  while(i<=d)
  {
    ith_bit=coeff(f,i);
    if(IsOne(ith_bit))
    {
      base=power(i); 
      num=num+base; 
    }   
    i++;
  }
  return num;
}

uint32_t power_uint32_t(long i)
{
  uint32_t num=1;
  for(long j=0;j<i;j++)
   num=num*2;
  return num;
}
//Conversion from GF2X to uint32_t

uint32_t GF2X_to_uint32_t(GF2X& f)
{
  uint32_t num=0;
  long d=deg(f);
  long i=0;
  uint32_t base=0;
  GF2 ith_bit;
  while(i<=d)
  {
    ith_bit=coeff(f,i);
    if(IsOne(ith_bit))
    {
      base=power_uint32_t(i); 
      num=num+base; 
    }   
    i++;
  }
  return num;
}

//Conversion of uint64_t to GF2X
GF2X uint64_t_to_GF2X(uint64_t num)
{
   GF2X f;
   f.SetLength(1000);
  
   long i=0;
   long ith_bit;
   long num_bits=NumBits(num);
        while(i<num_bits)
        {
          ith_bit=bit(num,i);
          SetCoeff(f,i,ith_bit);
          i++;
        }  
   f.normalize(); 
   return f;
 }

//Conversion of uint32_t to GF2X
GF2X uint32_t_to_GF2X(uint32_t num)
{
   GF2X f;
   f.SetLength(1000);
  
   long i=0;
   long ith_bit;
   long num_bits=NumBits(num);
        while(i<num_bits)
        {
          ith_bit=bit(num,i);
          SetCoeff(f,i,ith_bit);
          i++;
        }  
   f.normalize(); 
   return f;
 }

 void uint64_t_swap(uint64_t *x,uint64_t *y)
{
   uint64_t t;
   t=*x;
   *x=*y;
   *y=t;
} 
  
//Conversion of ZZ to polynomail over GF2
 GF2X set_GF2X(ZZ num)
{
  
  GF2X f;
  
  f.SetLength(1000); 
   
   long i=0;
   long ith_bit;

        while(!IsZero(num))
        {
          ith_bit=num % long(10);
          SetCoeff(f,i,ith_bit);
          i++;
          num=num/(long)10;
        }  
   f.normalize();
   return f;
 }

//helper function
ZZ func(char c,long length,long base,long i,ZZ val)
{
 ZZ temp1,temp3,temp4;
    long temp2;
    
    temp1=ZZ(c)-val;
    temp2=length-i;
    temp2=temp2-1;
    power(temp3,base,temp2);
    temp4=temp1*temp3;
    return temp4;
}

//hexadecimal to binary number for big integers
ZZ hex_to_binary(char hex[])
{
   ZZ num=hex_to_decimal(hex);   

   ZZ rem, i=ZZ(1);
   ZZ binary=ZZ(0);
  
    while (num!=0)
    {
        rem=num%2;
        num/=2;
        binary+=rem*i;
        i*=10;
    }
    return binary;
   
} 


//hexadecimal to decimal number for big integers
 ZZ hex_to_decimal(char *hex)
{
    char *hexString;
   long length = 0;
    long base = 16; // Base of Hexadecimal Number
    ZZ decimalNumber = ZZ(0);
    long i;
    
    
    // Find length of Hexadecimal Number
    for (hexString = hex; *hexString != '\0'; hexString++)
    {
        length++;
    }
    
    // Find Hexadecimal Number
    hexString = hex;
    
    for (i = 0; *hexString != '\0' && i < length; i++, hexString++)
    {
        // Compare *hexString with ASCII values
        if (*hexString >= 48 && *hexString <= 57)   // is *hexString Between 0-9
        {
             decimalNumber+=func(*hexString,length,base,i,ZZ(48));
        }
        else if ((*hexString >= 65 && *hexString <= 70))   // is *hexString Between A-F
        {
              decimalNumber+=func(*hexString,length,base,i,ZZ(55));
        }
        else if (*hexString >= 97 && *hexString <= 102)   // is *hexString Between a-f
        {
             decimalNumber+=func(*hexString,length,base,i,ZZ(87));
        }
      
    } 

     return decimalNumber;
     
   }


uint32_t power_uint32_t(int x,uint32_t i)
{ 
  uint32_t res=1;
  for(unsigned j=0;j<i;j++)
    res=res*x;
   return res;
}
 
//helper function
uint32_t func_uint32_t(char c,int length,int base,int i,int val)
{
 uint32_t temp1,temp3,temp4;
    int temp2;
    
    temp1=uint32_t(c)-val;
    temp2=length-i;
    temp2=temp2-1;
    temp3=power_uint32_t(base,temp2);
    temp4=temp1*temp3;
    return temp4;
}


//hexadecimal to uint32_t
uint32_t hex_to_uint32_t(char hex[])
{
    char *hexString;
    int length = 0;
    int base = 16; // Base of Hexadecimal Number
    uint32_t decimalNumber = 0;
    int i;
    
    
    // Find length of Hexadecimal Number
    for (hexString = hex; *hexString != '\0'; hexString++)
    {
        length++;
    }
    
    // Find Hexadecimal Number
    hexString = hex;
    
    for (i = 0; *hexString != '\0' && i < length; i++, hexString++)
    {
        // Compare *hexString with ASCII values
        if (*hexString >= 48 && *hexString <= 57)   // is *hexString Between 0-9
        {
             decimalNumber+=func_uint32_t(*hexString,length,base,i,(48));
        }
        else if ((*hexString >= 65 && *hexString <= 70))   // is *hexString Between A-F
        {
              decimalNumber+=func_uint32_t(*hexString,length,base,i,55);
        }
        else if (*hexString >= 97 && *hexString <= 102)   // is *hexString Between a-f
        {
             decimalNumber+=func_uint32_t(*hexString,length,base,i,87);
        }
      
    } 

     return decimalNumber;
     
   }