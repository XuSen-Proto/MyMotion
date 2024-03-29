function [dAiv1dr,dan,ddan] = dAiv1drchk(in1,in2)
%dAiv1drchk
%    [dAiv1dr,DAN,DDAN] = dAiv1drchk(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    08-Mar-2023 15:58:48

r1 = in1(1,:);
r2 = in1(2,:);
r3 = in1(3,:);
v11 = in2(1,:);
v12 = in2(2,:);
v13 = in2(3,:);
t2 = r1.^2;
t3 = r1.^3;
t4 = r2.^2;
t6 = r2.^3;
t7 = r3.^2;
t8 = r1.^5;
t10 = r3.^3;
t11 = r2.^5;
t13 = r3.^5;
t29 = r1.*r2.*r3.*v11.*4.0;
t30 = r1.*r2.*r3.*v12.*4.0;
t31 = r1.*r2.*r3.*v13.*4.0;
t5 = t2.^2;
t9 = t4.^2;
t12 = t7.^2;
t32 = t3.*v11.*2.0;
t33 = t6.*v12.*2.0;
t34 = t10.*v13.*2.0;
t35 = r1.*t4.*v11.*2.0;
t36 = r1.*t7.*v11.*2.0;
t37 = r2.*t2.*v12.*2.0;
t38 = r2.*t7.*v12.*2.0;
t39 = r3.*t2.*v13.*2.0;
t40 = r3.*t4.*v13.*2.0;
t41 = -t29;
t42 = -t30;
t43 = -t31;
t44 = r1.*r2.*t10.*v11;
t45 = r1.*r3.*t6.*v11;
t46 = r2.*r3.*t3.*v11;
t47 = r1.*r2.*t10.*v12;
t48 = r1.*r3.*t6.*v12;
t49 = r2.*r3.*t3.*v12;
t50 = r1.*r2.*t10.*v13;
t51 = r1.*r3.*t6.*v13;
t52 = r2.*r3.*t3.*v13;
t56 = t2.*t6.*v11;
t57 = t3.*t4.*v11;
t58 = t2.*t6.*v12;
t59 = t2.*t10.*v11;
t60 = t3.*t4.*v12;
t61 = t3.*t7.*v11;
t62 = t2.*t10.*v13;
t63 = t3.*t7.*v13;
t64 = t4.*t10.*v12;
t65 = t6.*t7.*v12;
t66 = t4.*t10.*v13;
t67 = t6.*t7.*v13;
t68 = t2.*t4.*v11.*2.0;
t69 = t2.*t4.*v12.*2.0;
t70 = t2.*t7.*v11.*2.0;
t71 = t2.*t4.*v13.*2.0;
t72 = t2.*t7.*v12.*2.0;
t73 = t4.*t7.*v11.*2.0;
t74 = t2.*t7.*v13.*2.0;
t75 = t4.*t7.*v12.*2.0;
t76 = t4.*t7.*v13.*2.0;
t77 = r1.*t4.*t7.*v11;
t78 = r2.*t2.*t7.*v11;
t79 = r3.*t2.*t4.*v11;
t80 = r1.*t4.*t7.*v12;
t81 = r2.*t2.*t7.*v12;
t82 = r3.*t2.*t4.*v12;
t83 = r1.*t4.*t7.*v13;
t84 = r2.*t2.*t7.*v13;
t85 = r3.*t2.*t4.*v13;
t86 = t2+t4+t7;
t14 = t5.*v11;
t15 = t5.*v12;
t16 = t9.*v11;
t17 = t5.*v13;
t18 = t9.*v12;
t19 = t12.*v11;
t20 = t9.*v13;
t21 = t12.*v12;
t22 = t12.*v13;
t53 = -t32;
t54 = -t33;
t55 = -t34;
t87 = -t77;
t88 = -t81;
t89 = -t85;
t90 = 1.0./t86.^2;
t91 = sqrt(t86);
t23 = r1.*t16;
t24 = r1.*t19;
t25 = r2.*t15;
t26 = r2.*t21;
t27 = r3.*t17;
t28 = r3.*t20;
t92 = t91.^3;
t93 = cos(t91);
t94 = sin(t91);
t95 = t93-1.0;
t96 = t14.*t93;
t97 = t15.*t93;
t98 = t16.*t93;
t99 = t17.*t93;
t100 = t18.*t93;
t101 = t19.*t93;
t102 = t20.*t93;
t103 = t21.*t93;
t104 = t22.*t93;
t105 = t32.*t93;
t106 = t33.*t93;
t107 = t34.*t93;
t109 = t29.*t93;
t110 = t30.*t93;
t111 = t31.*t93;
t112 = t35.*t93;
t113 = t36.*t93;
t114 = t37.*t93;
t115 = t38.*t93;
t116 = t39.*t93;
t117 = t40.*t93;
t118 = t68.*t93;
t119 = t69.*t93;
t120 = t70.*t93;
t121 = t71.*t93;
t122 = t72.*t93;
t123 = t73.*t93;
t124 = t74.*t93;
t125 = t75.*t93;
t126 = t76.*t93;
t127 = r1.*t92.*t94.*v11;
t128 = r2.*t92.*t94.*v12;
t129 = r3.*t92.*t94.*v13;
t130 = r1.*r2.*r3.*t91.*t94.*v11;
t131 = r1.*r2.*r3.*t91.*t94.*v12;
t132 = r1.*r2.*r3.*t91.*t94.*v13;
t133 = r1.*t4.*t91.*t94.*v11;
t134 = r1.*t7.*t91.*t94.*v11;
t135 = r2.*t2.*t91.*t94.*v12;
t136 = r2.*t7.*t91.*t94.*v12;
t137 = r3.*t2.*t91.*t94.*v13;
t138 = r3.*t4.*t91.*t94.*v13;
t108 = 1.0./t95;
mt1 = [(t90.*t108.*(t23+t24-t25-t27+t37-t38+t39-t40+t54+t55+t57-t58+t61-t62+t88+t89+t106+t107+t115+t117+t128+t129+t133+t134-t135-t137+t7.*t35-r1.*t4.*v11.*4.0-r1.*t7.*v11.*4.0+r1.*t4.*t93.*v11.*4.0-r2.*t2.*t93.*v12.*2.0-r3.*t2.*t93.*v13.*2.0+r1.*t7.*t93.*v11.*4.0))./2.0];
mt2 = [t90.*t108.*(-t17-t20-t22+t43+t50+t51+t52+t56-t60-t71-t74-t76+t78-t80+t99+t102+t104+t111+t121+t124+t126+t132+r2.*t14-r1.*t21+t6.*v11.*2.0-t8.*v12-r2.*t2.*v11.*2.0-r1.*t4.*v12.*4.0+r2.*t7.*v11.*2.0-t3.*t7.*v12.*2.0-t6.*t93.*v11.*2.0+r2.*t2.*t93.*v11.*2.0+r1.*t4.*t93.*v12.*4.0-r2.*t7.*t93.*v11.*2.0+r1.*t92.*t94.*v12.*2.0-r2.*t92.*t94.*v11-t3.*t91.*t94.*v12+r2.*t2.*t91.*t94.*v11-r1.*t7.*t91.*t94.*v12).*(-1.0./2.0)];
mt3 = [t90.*t108.*(t15+t18+t21+t42+t47+t48+t49+t59-t63+t69+t72+t75+t79-t83-t97-t100-t103+t110+t131+r3.*t14-r1.*t20-t8.*v13+t10.*v11.*2.0-r3.*t2.*v11.*2.0+r3.*t4.*v11.*2.0-r1.*t7.*v13.*4.0-t3.*t4.*v13.*2.0-t10.*t93.*v11.*2.0+r3.*t2.*t93.*v11.*2.0-r3.*t4.*t93.*v11.*2.0+r1.*t7.*t93.*v13.*4.0+r1.*t92.*t94.*v13.*2.0-r3.*t92.*t94.*v11-t2.*t4.*t93.*v12.*2.0-t2.*t7.*t93.*v12.*2.0-t4.*t7.*t93.*v12.*2.0-t3.*t91.*t94.*v13+r3.*t2.*t91.*t94.*v11-r1.*t4.*t91.*t94.*v13).*(-1.0./2.0)];
mt4 = [t90.*t108.*(t17+t20+t22+t43+t50+t51+t52-t56+t60+t71+t74+t76-t78+t80-t99-t102-t104+t111+t132+r1.*t18-r2.*t19+t3.*v12.*2.0-t11.*v11-r2.*t2.*v11.*4.0-r1.*t4.*v12.*2.0+r1.*t7.*v12.*2.0-t6.*t7.*v11.*2.0-t3.*t93.*v12.*2.0+r2.*t2.*t93.*v11.*4.0+r1.*t4.*t93.*v12.*2.0-r1.*t7.*t93.*v12.*2.0-r1.*t92.*t94.*v12+r2.*t92.*t94.*v11.*2.0-t2.*t4.*t93.*v13.*2.0-t2.*t7.*t93.*v13.*2.0-t4.*t7.*t93.*v13.*2.0-t6.*t91.*t94.*v11+r1.*t4.*t91.*t94.*v12-r2.*t7.*t91.*t94.*v11).*(-1.0./2.0)];
mt5 = [(t90.*t108.*(-t23+t25+t26-t28+t35-t36-t39+t40+t53+t55-t57+t58+t65-t66+t87+t89+t105+t107+t113+t116+t127+t129-t133+t135+t136-t138+t7.*t37-r2.*t2.*v12.*4.0-r2.*t7.*v12.*4.0-r1.*t4.*t93.*v11.*2.0+r2.*t2.*t93.*v12.*4.0-r3.*t4.*t93.*v13.*2.0+r2.*t7.*t93.*v12.*4.0))./2.0];
mt6 = [t90.*t108.*(-t14-t16-t19+t41+t44+t45+t46+t64-t67-t68-t70-t73+t82-t84+t96+t98+t101+t109+t118+t120+t123+t130-r2.*t17+r3.*t18+t10.*v12.*2.0-t11.*v13+r3.*t2.*v12.*2.0-r3.*t4.*v12.*2.0-r2.*t7.*v13.*4.0-t2.*t6.*v13.*2.0-t10.*t93.*v12.*2.0-r3.*t2.*t93.*v12.*2.0+r3.*t4.*t93.*v12.*2.0+r2.*t7.*t93.*v13.*4.0+r2.*t92.*t94.*v13.*2.0-r3.*t92.*t94.*v12-t6.*t91.*t94.*v13-r2.*t2.*t91.*t94.*v13+r3.*t4.*t91.*t94.*v12).*(-1.0./2.0)];
mt7 = [t90.*t108.*(-t15-t18-t21+t42+t47+t48+t49-t59+t63-t69-t72-t75-t79+t83+t97+t100+t103+t110+t119+t122+t125+t131-r3.*t16+r1.*t22+t3.*v13.*2.0-t13.*v11-r3.*t2.*v11.*4.0+r1.*t4.*v13.*2.0-r1.*t7.*v13.*2.0-t4.*t10.*v11.*2.0-t3.*t93.*v13.*2.0+r3.*t2.*t93.*v11.*4.0-r1.*t4.*t93.*v13.*2.0+r1.*t7.*t93.*v13.*2.0-r1.*t92.*t94.*v13+r3.*t92.*t94.*v11.*2.0-t10.*t91.*t94.*v11-r3.*t4.*t91.*t94.*v11+r1.*t7.*t91.*t94.*v13).*(-1.0./2.0)];
mt8 = [t90.*t108.*(t14+t16+t19+t41+t44+t45+t46-t64+t67+t68+t70+t73-t82+t84-t96-t98-t101+t109+t130-r3.*t15+r2.*t22+t6.*v13.*2.0-t13.*v12+r2.*t2.*v13.*2.0-r3.*t4.*v12.*4.0-r2.*t7.*v13.*2.0-t2.*t10.*v12.*2.0-t6.*t93.*v13.*2.0-r2.*t2.*t93.*v13.*2.0+r3.*t4.*t93.*v12.*4.0+r2.*t7.*t93.*v13.*2.0-r2.*t92.*t94.*v13+r3.*t92.*t94.*v12.*2.0-t2.*t4.*t93.*v11.*2.0-t2.*t7.*t93.*v11.*2.0-t4.*t7.*t93.*v11.*2.0-t10.*t91.*t94.*v12-r3.*t2.*t91.*t94.*v12+r2.*t7.*t91.*t94.*v13).*(-1.0./2.0)];
mt9 = [(t90.*t108.*(-t24-t26+t27+t28-t35+t36-t37+t38+t53+t54-t61+t62-t65+t66+t87+t88+t105+t106+t112+t114+t127+t128-t134-t136+t137+t138+t4.*t39-r3.*t2.*v13.*4.0-r3.*t4.*v13.*4.0+r3.*t2.*t93.*v13.*4.0-r1.*t7.*t93.*v11.*2.0+r3.*t4.*t93.*v13.*4.0-r2.*t7.*t93.*v12.*2.0))./2.0];
dAiv1dr = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9],3,3);
if nargout > 1
    dan = [0.0,0.0,t108.*(t86+t95.*4.0+t91.*t94).*(-1.0./2.0)];
end
if nargout > 2
    ddan = [0.0,0.0,t108.^2.*(t86.*3.0+t93.*1.6e+1+t95.*1.6e+1-t86.*t93.*3.0+t91.*t94.*3.0+t92.*t94-t93.^2.*1.6e+1-t91.*t93.*t94.*3.0).*(-1.0./2.0)];
end
