function T1T2 = SE3x(in1,in2)
%SE3x
%    T1T2 = SE3x(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    2024-01-16 21:47:14

T1_1_1 = in1(1);
T1_1_2 = in1(5);
T1_1_3 = in1(9);
T1_1_4 = in1(13);
T1_2_1 = in1(2);
T1_2_2 = in1(6);
T1_2_3 = in1(10);
T1_2_4 = in1(14);
T1_3_1 = in1(3);
T1_3_2 = in1(7);
T1_3_3 = in1(11);
T1_3_4 = in1(15);
T2_1_1 = in2(1);
T2_1_2 = in2(5);
T2_1_3 = in2(9);
T2_1_4 = in2(13);
T2_2_1 = in2(2);
T2_2_2 = in2(6);
T2_2_3 = in2(10);
T2_2_4 = in2(14);
T2_3_1 = in2(3);
T2_3_2 = in2(7);
T2_3_3 = in2(11);
T2_3_4 = in2(15);
T1T2 = reshape([T1_1_1.*T2_1_1+T1_1_2.*T2_2_1+T1_1_3.*T2_3_1,T1_2_1.*T2_1_1+T1_2_2.*T2_2_1+T1_2_3.*T2_3_1,T1_3_1.*T2_1_1+T1_3_2.*T2_2_1+T1_3_3.*T2_3_1,0.0,T1_1_1.*T2_1_2+T1_1_2.*T2_2_2+T1_1_3.*T2_3_2,T1_2_1.*T2_1_2+T1_2_2.*T2_2_2+T1_2_3.*T2_3_2,T1_3_1.*T2_1_2+T1_3_2.*T2_2_2+T1_3_3.*T2_3_2,0.0,T1_1_1.*T2_1_3+T1_1_2.*T2_2_3+T1_1_3.*T2_3_3,T1_2_1.*T2_1_3+T1_2_2.*T2_2_3+T1_2_3.*T2_3_3,T1_3_1.*T2_1_3+T1_3_2.*T2_2_3+T1_3_3.*T2_3_3,0.0,T1_1_4+T1_1_1.*T2_1_4+T1_1_2.*T2_2_4+T1_1_3.*T2_3_4,T1_2_4+T1_2_1.*T2_1_4+T1_2_2.*T2_2_4+T1_2_3.*T2_3_4,T1_3_4+T1_3_1.*T2_1_4+T1_3_2.*T2_2_4+T1_3_3.*T2_3_4,1.0],[4,4]);
end
