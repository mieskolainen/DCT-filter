#include "dct.h"


/*2D forward DCT transform*/
void fdct(const double *inblock, double *outblock)
{
   
     double s[8],t[8],r[8];
     const double *p;
     int j;
     /*Horizontal direction*/
     p=inblock;
     for(j=0;j<64;j+=8)
     {
	
	      /*First stage*/
	      s[0]=(*(p)+*(p+7));
	      s[1]=(*(p+1)+*(p+6));
	      s[2]=(*(p+2)+*(p+5));
	      s[3]=(*(p+3)+*(p+4));
	      s[4]=(*(p+3)-*(p+4));
	      s[5]=(*(p+2)-*(p+5));
	      s[6]=(*(p+1)-*(p+6));
	      s[7]=(*(p)-*(p+7));
	      /*Second stage*/
	      t[0]=s[0]+s[3];
	      t[1]=s[1]+s[2];
	      t[2]=s[1]-s[2];
	      t[3]=s[0]-s[3];
	      t[5]=(s[6]-s[5])*W3;
	      t[6]=(s[6]+s[5])*W3;
	      /*Third stage*/
	      r[4]=s[4]+t[5];
	      r[5]=s[4]-t[5];
	      r[6]=s[7]-t[6];
	      r[7]=s[7]+t[6];
	      /*Fourth stage*/
	      outblock[j]=(t[0]+t[1])*W3;
	      outblock[4+j]=(t[0]-t[1])*W3;
	      outblock[1+j]=(r[4]*W5+r[7]*W4);
	      outblock[7+j]=(r[7]*W5-r[4]*W4);
	      outblock[3+j]=(r[6]*W6-r[5]*W7);
	      outblock[5+j]=(r[5]*W6+r[6]*W7);
	      outblock[2+j]=(t[2]*W2+t[3]*W1);
	      outblock[6+j]=(t[3]*W2-t[2]*W1);
	      p +=8;
     }
   
         /*Vertical direction*/
     for(j=0;j<8;j++)
     {
	
	      /*First stage*/
	      s[0]=outblock[j]+outblock[j+56];
	      s[1]=outblock[j+8]+outblock[j+48];
	      s[2]=outblock[j+16]+outblock[j+40];
	      s[3]=outblock[j+24]+outblock[j+32];
	      s[4]=outblock[j+24]-outblock[j+32];
	      s[5]=outblock[j+16]-outblock[j+40];
	      s[6]=outblock[j+8]-outblock[j+48];
	      s[7]=outblock[j]-outblock[j+56];
	      /*Second stage*/
	      t[0]=s[0]+s[3];
	      t[1]=s[1]+s[2];
	      t[2]=s[1]-s[2];
	      t[3]=s[0]-s[3];
	      t[5]=(s[6]-s[5])*W3;
	      t[6]=(s[6]+s[5])*W3;
	      /*Third stage*/
	      r[4]=s[4]+t[5];
	      r[5]=s[4]-t[5];
	      r[6]=s[7]-t[6];
	      r[7]=s[7]+t[6];
	      /*Fourth stage, transform coefficients scaled by 2/N*/
	      outblock[j]=(t[0]+t[1])*W3/4.0;
	      outblock[32+j]=(t[0]-t[1])*W3/4.0;
	      outblock[8+j]=(r[4]*W5+r[7]*W4)/4.0;
	      outblock[56+j]=(r[7]*W5-r[4]*W4)/4.0;
	      outblock[24+j]=(r[6]*W6-r[5]*W7)/4.0;
	      outblock[40+j]=(r[5]*W6+r[6]*W7)/4.0;
	      outblock[16+j]=(t[2]*W2+t[3]*W1)/4.0;
	      outblock[48+j]=(t[3]*W2-t[2]*W1)/4.0;
     }
   
}


/*2D inverse DCT transfrom*/
void idct(const double *inblock, double *outblock) 
{
   
     double s[8],t[8],r[8];
     const double *p;
     int j;
     /*Horizontal direction*/
     p=inblock;
     for(j=0;j<64;j+=8)
     {
	
	      /*First stage*/
	      r[4]=*(p+1)*W5-*(p+7)*W4;
	      r[5]=*(p+5)*W6-*(p+3)*W7;
	      r[6]=*(p+3)*W6+*(p+5)*W7;
	      r[7]=*(p+7)*W5+*(p+1)*W4;
	      
	      /*Second stage*/
	      t[0]=(*(p)+*(p+4))*W3;
	      t[1]=(*(p)-*(p+4))*W3;
	      t[2]=*(p+2)*W2-*(p+6)*W1;
	      t[3]=*(p+2)*W1+*(p+6)*W2;
	      t[4]=(r[4]+r[5]);
	      t[5]=(r[4]-r[5]);
	      t[6]=(r[7]-r[6]);
	      t[7]=(r[6]+r[7]);
	      
	      /*Third stage*/
	      s[0]=(t[0]+t[3]);
	      s[1]=(t[1]+t[2]);
	      s[2]=(t[1]-t[2]);
	      s[3]=(t[0]-t[3]);
	      s[5]=(t[6]-t[5])*W3;
	      s[6]=(t[5]+t[6])*W3;
	      
	      /*Fourth stage*/
	      outblock[j]=(s[0]+t[7]);
	      outblock[j+1]=(s[1]+s[6]);
	      outblock[j+2]=(s[2]+s[5]);
	      outblock[j+3]=(s[3]+t[4]);
	      outblock[j+4]=(s[3]-t[4]);
	      outblock[j+5]=(s[2]-s[5]);
	      outblock[j+6]=(s[1]-s[6]);
	      outblock[j+7]=(s[0]-t[7]);
	      p +=8;
     }
   
     /*Vertical direction*/
     for(j=0;j<8;j++)
     {
	
	      /*First stage*/
	      r[4]=outblock[j+8]*W5-outblock[j+56]*W4;
	      r[5]=outblock[j+40]*W6-outblock[j+24]*W7;
	      r[6]=outblock[j+24]*W6+outblock[j+40]*W7;
	      r[7]=outblock[j+56]*W5+outblock[j+8]*W4;
	
	      /*Second stage*/
	      t[0]=(outblock[j]+outblock[j+32])*W3;
	      t[1]=(outblock[j]-outblock[j+32])*W3;
	      t[2]=outblock[j+16]*W2-outblock[j+48]*W1;
	      t[3]=outblock[j+16]*W1+outblock[j+48]*W2;
	      t[4]=(r[4]+r[5]);
	      t[5]=(r[4]-r[5]);
	      t[6]=(r[7]-r[6]);
	      t[7]=(r[6]+r[7]);
	
	      /*Third stage*/
	      s[0]=(t[0]+t[3]);
	      s[1]=(t[1]+t[2]);
	      s[2]=(t[1]-t[2]);
	      s[3]=(t[0]-t[3]);
	      s[5]=(t[6]-t[5])*W3;
	      s[6]=(t[5]+t[6])*W3;
	
	      /*Fourth stage, pixel values scaled by 2/N*/
	      outblock[j]=(s[0]+t[7])/4.0;
	      outblock[j+8]=(s[1]+s[6])/4.0;
	      outblock[j+16]=(s[2]+s[5])/4.0;
	      outblock[j+24]=(s[3]+t[4])/4.0;
	      outblock[j+32]=(s[3]-t[4])/4.0;
	      outblock[j+40]=(s[2]-s[5])/4.0;
	      outblock[j+48]=(s[1]-s[6])/4.0;
	      outblock[j+56]=(s[0]-t[7])/4.0;
     }
   
}
