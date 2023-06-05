function [Center,Boundry1,Boundry2,Boundry3,Boundry4,Boundry5,Boundry6]=Creat_Periphery_Mask(Mask,degree)
   
[Periphery_Inside1 Periphery_Outside1 Core_1]=FindPeriphery(Mask, degree);
MaskPWS_Periphery_Inside1=Periphery_Inside1;
MaskPWS_Periphery_Inside1(MaskPWS_Periphery_Inside1>0)=1;

[Periphery_Inside2 Periphery_Outside2 Core_2]=FindPeriphery(Mask, 2*degree);
MaskPWS_Periphery_Inside2=Periphery_Inside2-Periphery_Inside1;
MaskPWS_Periphery_Inside2(MaskPWS_Periphery_Inside2>0)=1;
MaskPWS_Periphery_Inside2(MaskPWS_Periphery_Inside2<0)=0;

[Periphery_Inside3 Periphery_Outside3 Core_3]=FindPeriphery(Mask, 3*degree);
MaskPWS_Periphery_Inside3=Periphery_Inside3-Periphery_Inside2;
MaskPWS_Periphery_Inside3(MaskPWS_Periphery_Inside3>0)=1;
MaskPWS_Periphery_Inside3(MaskPWS_Periphery_Inside3<0)=0;


[Periphery_Inside4 Periphery_Outside4 Core_4]=FindPeriphery(Mask, 4*degree);
MaskPWS_Periphery_Inside4=Periphery_Inside4-Periphery_Inside3;
MaskPWS_Periphery_Inside4(MaskPWS_Periphery_Inside4>0)=1;
MaskPWS_Periphery_Inside4(MaskPWS_Periphery_Inside4<0)=0;


[Periphery_Inside5 Periphery_Outside5 Core_5]=FindPeriphery(Mask, 5*degree);
MaskPWS_Periphery_Inside5=Periphery_Inside5-Periphery_Inside4;
MaskPWS_Periphery_Inside5(MaskPWS_Periphery_Inside5>0)=1;
MaskPWS_Periphery_Inside5(MaskPWS_Periphery_Inside5<0)=0;


[Periphery_Inside6 Periphery_Outside6 Core_6]=FindPeriphery(Mask, 6*degree);
MaskPWS_Periphery_Inside6=Periphery_Inside6-Periphery_Inside5;
MaskPWS_Periphery_Inside6(MaskPWS_Periphery_Inside6>0)=1;
MaskPWS_Periphery_Inside6(MaskPWS_Periphery_Inside6<0)=0;


Boundry1=MaskPWS_Periphery_Inside1;
Boundry2=MaskPWS_Periphery_Inside2;
Boundry3=MaskPWS_Periphery_Inside3;
Boundry4=MaskPWS_Periphery_Inside4;
Boundry5=MaskPWS_Periphery_Inside5;
Boundry6=MaskPWS_Periphery_Inside6;


[Periphery_Inside8 Periphery_Outside8 Core_8]=FindPeriphery(Mask, 8*degree);
Center=Core_8;



