function [dAv1drTv2_dr,dAv1drTv2_dv1,C] = dwbdso3chk(in1,in2,in3)
%DWBDSO3CHK
%    [dAv1drTv2_dr,dAv1drTv2_dv1,C] = DWBDSO3CHK(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    08-Mar-2023 19:19:56

r1 = in1(1,:);
r2 = in1(2,:);
r3 = in1(3,:);
v11 = in2(1,:);
v12 = in2(2,:);
v13 = in2(3,:);
v21 = in3(1,:);
v22 = in3(2,:);
v23 = in3(3,:);
t2 = r1.^2;
t3 = r1.^3;
t4 = r2.^2;
t6 = r2.^3;
t7 = r3.^2;
t9 = r3.^3;
t5 = t2.^2;
t8 = t4.^2;
t10 = t7.^2;
t11 = t2+t4+t7;
t12 = t11.^2;
t13 = 1.0./t11;
t15 = sqrt(t11);
t14 = 1.0./t12;
t16 = t15.^3;
t17 = t15.^5;
t18 = 1.0./t15;
t22 = cos(t15);
t23 = sin(t15);
t27 = t15.*2.0;
t35 = -t15;
t45 = r1.*r2.*t15.*-2.0;
t46 = r1.*r3.*t15.*-2.0;
t47 = r2.*r3.*t15.*-2.0;
t96 = r1.*r2.*r3.*t15.*v11.*v21.*8.0;
t97 = r1.*r2.*r3.*t15.*v12.*v22.*8.0;
t98 = r1.*r2.*r3.*t15.*v13.*v23.*8.0;
t99 = r1.*t6.*t15.*v11.*v22.*8.0;
t100 = r2.*t2.*t15.*v13.*v21.*8.0;
t101 = r3.*t2.*t15.*v12.*v21.*8.0;
t102 = r2.*t3.*t15.*v12.*v21.*8.0;
t103 = r1.*t4.*t15.*v13.*v22.*8.0;
t104 = r3.*t4.*t15.*v11.*v22.*8.0;
t105 = r1.*t7.*t15.*v12.*v23.*8.0;
t106 = r1.*t9.*t15.*v11.*v23.*8.0;
t107 = r2.*t7.*t15.*v11.*v23.*8.0;
t108 = r3.*t3.*t15.*v13.*v21.*8.0;
t109 = r2.*t9.*t15.*v12.*v23.*8.0;
t110 = r3.*t6.*t15.*v13.*v22.*8.0;
t111 = r2.*r3.*t2.*t15.*v11.*v21.*8.0;
t112 = r1.*r2.*t7.*t15.*v11.*v22.*8.0;
t113 = r1.*r2.*t7.*t15.*v12.*v21.*8.0;
t117 = r1.*r3.*t4.*t15.*v11.*v23.*8.0;
t118 = r1.*r3.*t4.*t15.*v12.*v22.*8.0;
t119 = r1.*r3.*t4.*t15.*v13.*v21.*8.0;
t123 = r2.*r3.*t2.*t15.*v12.*v23.*8.0;
t124 = r2.*r3.*t2.*t15.*v13.*v22.*8.0;
t125 = r1.*r2.*t7.*t15.*v13.*v23.*8.0;
t132 = t2.*t4.*t15.*v11.*v21.*8.0;
t133 = t2.*t7.*t15.*v11.*v21.*8.0;
t134 = t2.*t4.*t15.*v12.*v22.*8.0;
t135 = t4.*t7.*t15.*v12.*v22.*8.0;
t136 = t2.*t7.*t15.*v13.*v23.*8.0;
t137 = t4.*t7.*t15.*v13.*v23.*8.0;
t139 = r1.*r3.*t4.*t15.*v11.*v22.*-8.0;
t140 = r2.*r3.*t2.*t15.*v12.*v21.*-8.0;
t141 = r1.*r2.*t7.*t15.*v11.*v23.*-8.0;
t143 = r2.*r3.*t2.*t15.*v13.*v21.*-8.0;
t144 = r1.*r2.*t7.*t15.*v12.*v23.*-8.0;
t145 = r1.*r3.*t4.*t15.*v13.*v22.*-8.0;
t19 = 1.0./t16;
t20 = 1.0./t17;
t21 = t18.^7;
t24 = r1.*t16;
t25 = r2.*t16;
t26 = r3.*t16;
t28 = r1.*t18;
t29 = r2.*t18;
t30 = r3.*t18;
t31 = t23.*3.0;
t32 = t17.*v11.*v21;
t33 = t17.*v12.*v22;
t34 = t17.*v13.*v23;
t36 = t22.*2.0;
t38 = r1.*r2.*t27;
t39 = r1.*r3.*t27;
t40 = r2.*r3.*t27;
t44 = t22-1.0;
t63 = t3.*t23.*v21;
t64 = t6.*t23.*v22;
t65 = t9.*t23.*v23;
t75 = t3.*t27.*v21;
t76 = t6.*t27.*v22;
t77 = t9.*t27.*v23;
t78 = r1.*t4.*t27.*v21;
t79 = r1.*t7.*t27.*v21;
t80 = r2.*t2.*t27.*v22;
t81 = r2.*t7.*t27.*v22;
t82 = r3.*t2.*t27.*v23;
t83 = r3.*t4.*t27.*v23;
t84 = t2.*t16.*v11.*v21.*2.0;
t85 = t4.*t16.*v11.*v21.*2.0;
t86 = t2.*t16.*v12.*v22.*2.0;
t87 = t7.*t16.*v11.*v21.*2.0;
t88 = t4.*t16.*v12.*v22.*2.0;
t89 = t2.*t16.*v13.*v23.*2.0;
t90 = t7.*t16.*v12.*v22.*2.0;
t91 = t4.*t16.*v13.*v23.*2.0;
t92 = t7.*t16.*v13.*v23.*2.0;
t114 = r1.*t104;
t115 = r2.*t101;
t116 = r1.*t107;
t120 = r3.*t100;
t121 = r2.*t105;
t122 = r3.*t103;
t138 = -t111;
t142 = -t118;
t146 = -t125;
t147 = r2.*r3.*t2.*t23.*v11.*v21.*1.5e+1;
t148 = r1.*r2.*t7.*t23.*v11.*v22.*1.5e+1;
t149 = r1.*r2.*t7.*t23.*v12.*v21.*1.5e+1;
t150 = r1.*r3.*t4.*t23.*v11.*v22.*1.5e+1;
t151 = r2.*r3.*t2.*t23.*v12.*v21.*1.5e+1;
t152 = r1.*r2.*t7.*t23.*v11.*v23.*1.5e+1;
t153 = r1.*r3.*t4.*t23.*v11.*v23.*1.5e+1;
t154 = r1.*r3.*t4.*t23.*v12.*v22.*1.5e+1;
t155 = r1.*r3.*t4.*t23.*v13.*v21.*1.5e+1;
t156 = r2.*r3.*t2.*t23.*v13.*v21.*1.5e+1;
t157 = r1.*r2.*t7.*t23.*v12.*v23.*1.5e+1;
t158 = r1.*r3.*t4.*t23.*v13.*v22.*1.5e+1;
t159 = r2.*r3.*t2.*t23.*v12.*v23.*1.5e+1;
t160 = r2.*r3.*t2.*t23.*v13.*v22.*1.5e+1;
t161 = r1.*r2.*t7.*t23.*v13.*v23.*1.5e+1;
t162 = t2.*t4.*t23.*v11.*v21.*1.5e+1;
t163 = t2.*t7.*t23.*v11.*v21.*1.5e+1;
t164 = t2.*t4.*t23.*v12.*v22.*1.5e+1;
t165 = t4.*t7.*t23.*v12.*v22.*1.5e+1;
t166 = t2.*t7.*t23.*v13.*v23.*1.5e+1;
t167 = t4.*t7.*t23.*v13.*v23.*1.5e+1;
t168 = r1.*t11.*t23;
t169 = r2.*t11.*t23;
t170 = r3.*t11.*t23;
t174 = t12.*t23.*v11.*v21;
t175 = t12.*t23.*v12.*v22;
t176 = t12.*t23.*v13.*v23;
t177 = t15.*t22;
t178 = t15.*t23;
t180 = r1.*t12.*t23.*v12.*v23;
t181 = r1.*t12.*t23.*v13.*v22;
t182 = r2.*t12.*t23.*v11.*v23;
t184 = r2.*t12.*t23.*v13.*v21;
t185 = r3.*t12.*t23.*v11.*v22;
t186 = r3.*t12.*t23.*v12.*v21;
t255 = t2.*t11.*t23.*v11.*v21.*-3.0;
t256 = t4.*t11.*t23.*v12.*v22.*-3.0;
t257 = t7.*t11.*t23.*v13.*v23.*-3.0;
t258 = t23+t35;
t262 = t2.*t4.*t11.*t23.*v11.*v21;
t263 = t2.*t7.*t11.*t23.*v11.*v21;
t264 = t2.*t4.*t11.*t23.*v12.*v22;
t265 = t4.*t7.*t11.*t23.*v12.*v22;
t266 = t2.*t7.*t11.*t23.*v13.*v23;
t267 = t4.*t7.*t11.*t23.*v13.*v23;
t314 = t2.*t16.*t22.*v11.*v21;
t315 = t4.*t16.*t22.*v11.*v21;
t316 = t2.*t16.*t22.*v12.*v22;
t317 = t7.*t16.*t22.*v11.*v21;
t318 = t4.*t16.*t22.*v12.*v22;
t319 = t2.*t16.*t22.*v13.*v23;
t320 = t7.*t16.*t22.*v12.*v22;
t321 = t4.*t16.*t22.*v13.*v23;
t322 = t7.*t16.*t22.*v13.*v23;
t338 = t22.*t96;
t339 = t22.*t97;
t340 = t22.*t98;
t341 = t22.*t100;
t342 = t22.*t101;
t343 = t22.*t103;
t344 = t22.*t104;
t345 = t22.*t105;
t346 = t22.*t107;
t37 = -t31;
t41 = t24.*v21.*2.0;
t42 = t25.*v22.*2.0;
t43 = t26.*v23.*2.0;
t51 = -t32;
t52 = -t33;
t53 = -t34;
t55 = t24.*v12.*v23.*2.0;
t56 = t24.*v13.*v22.*2.0;
t57 = t25.*v11.*v23.*2.0;
t59 = t25.*v13.*v21.*2.0;
t60 = t26.*v11.*v22.*2.0;
t61 = t26.*v12.*v21.*2.0;
t66 = r2.*t24.*v11.*v23.*2.0;
t67 = r3.*t24.*v11.*v22.*2.0;
t68 = r3.*t25.*v11.*v21.*2.0;
t70 = r3.*t24.*v12.*v22.*2.0;
t71 = r3.*t25.*v12.*v21.*2.0;
t72 = r2.*t24.*v13.*v23.*2.0;
t93 = -t63;
t94 = -t64;
t95 = -t65;
t126 = r1.*t64.*v11.*1.5e+1;
t127 = r2.*t63.*v12.*1.5e+1;
t128 = r1.*t65.*v11.*1.5e+1;
t129 = r3.*t63.*v13.*1.5e+1;
t130 = r2.*t65.*v12.*1.5e+1;
t131 = r3.*t64.*v13.*1.5e+1;
t171 = r2.*t168;
t172 = r3.*t168;
t173 = r3.*t169;
t179 = r1.*t174;
t183 = r2.*t175;
t187 = r3.*t176;
t188 = r1.*r2.*t11.*t31.*v11.*v23;
t189 = r1.*r3.*t11.*t31.*v11.*v22;
t190 = r2.*r3.*t11.*t31.*v11.*v21;
t191 = r1.*r2.*t11.*t31.*v12.*v23;
t192 = r1.*r3.*t11.*t31.*v12.*v22;
t193 = r2.*r3.*t11.*t31.*v12.*v21;
t194 = r1.*r2.*t11.*t31.*v13.*v23;
t195 = r1.*r3.*t11.*t31.*v13.*v22;
t196 = r2.*r3.*t11.*t31.*v13.*v21;
t197 = r1.*t11.*t64.*v11;
t198 = r2.*t11.*t63.*v12;
t199 = r1.*t11.*t65.*v11;
t200 = r3.*t11.*t63.*v13;
t201 = r2.*t11.*t65.*v12;
t202 = r3.*t11.*t64.*v13;
t206 = t13.*t44;
t207 = t2.*t11.*t31.*v11.*v21;
t208 = t4.*t11.*t31.*v11.*v21;
t209 = t2.*t11.*t31.*v12.*v22;
t210 = t7.*t11.*t31.*v11.*v21;
t211 = t4.*t11.*t31.*v12.*v22;
t212 = t2.*t11.*t31.*v13.*v23;
t213 = t7.*t11.*t31.*v12.*v22;
t214 = t4.*t11.*t31.*v13.*v23;
t215 = t7.*t11.*t31.*v13.*v23;
t216 = t22.*t28;
t217 = t22.*t29;
t218 = t22.*t30;
t228 = t2.*t169.*v13.*v21.*5.0;
t229 = t2.*t170.*v12.*v21.*5.0;
t230 = t4.*t168.*v13.*v22.*5.0;
t231 = t4.*t170.*v11.*v22.*5.0;
t232 = t7.*t168.*v12.*v23.*5.0;
t233 = t7.*t169.*v11.*v23.*5.0;
t234 = t22.*t38;
t235 = t22.*t39;
t236 = t22.*t40;
t237 = r1.*r2.*t19.*t23;
t238 = r1.*r3.*t19.*t23;
t239 = r2.*r3.*t19.*t23;
t259 = t2.*t19.*t23;
t260 = t4.*t19.*t23;
t261 = t7.*t19.*t23;
t268 = t3.*t177.*v21;
t269 = t6.*t177.*v22;
t270 = t9.*t177.*v23;
t271 = t24.*t36.*v11.*v21;
t272 = t24.*t36.*v12.*v23;
t273 = t24.*t36.*v13.*v22;
t274 = t25.*t36.*v11.*v23;
t275 = t25.*t36.*v12.*v22;
t276 = t25.*t36.*v13.*v21;
t277 = t26.*t36.*v11.*v22;
t278 = t26.*t36.*v12.*v21;
t279 = t26.*t36.*v13.*v23;
t284 = r2.*t22.*t24.*v11.*v23;
t285 = r3.*t22.*t24.*v11.*v22;
t286 = r3.*t22.*t25.*v11.*v21;
t287 = r2.*t22.*t24.*v12.*v23;
t288 = r3.*t22.*t24.*v12.*v22;
t289 = r3.*t22.*t25.*v12.*v21;
t290 = r2.*t22.*t24.*v13.*v23;
t291 = r3.*t22.*t24.*v13.*v22;
t292 = r3.*t22.*t25.*v13.*v21;
t293 = r2.*r3.*t22.*t24.*v11.*v21;
t296 = r1.*r2.*t14.*t44.*2.0;
t297 = r1.*r3.*t14.*t44.*2.0;
t298 = r2.*r3.*t14.*t44.*2.0;
t308 = r1.*t4.*t177.*v21;
t309 = r1.*t7.*t177.*v21;
t310 = r2.*t2.*t177.*v22;
t311 = r2.*t7.*t177.*v22;
t312 = r3.*t2.*t177.*v23;
t313 = r3.*t4.*t177.*v23;
t323 = t2.*t22.*t25.*v13.*v21;
t324 = t2.*t22.*t26.*v12.*v21;
t325 = t4.*t22.*t24.*v13.*v22;
t326 = t4.*t22.*t26.*v11.*v22;
t327 = t7.*t22.*t24.*v12.*v23;
t328 = t7.*t22.*t25.*v11.*v23;
t329 = t2.*t14.*t44.*2.0;
t330 = t4.*t14.*t44.*2.0;
t331 = t7.*t14.*t44.*2.0;
t347 = r2.*r3.*t2.*t177.*v11.*v21.*7.0;
t350 = r1.*r3.*t4.*t177.*v11.*v22.*7.0;
t351 = r2.*r3.*t2.*t177.*v12.*v21.*7.0;
t352 = r1.*r2.*t7.*t177.*v11.*v23.*7.0;
t354 = r1.*r3.*t4.*t177.*v12.*v22.*7.0;
t356 = r2.*r3.*t2.*t177.*v13.*v21.*7.0;
t357 = r1.*r2.*t7.*t177.*v12.*v23.*7.0;
t358 = r1.*r3.*t4.*t177.*v13.*v22.*7.0;
t361 = r1.*r2.*t7.*t177.*v13.*v23.*7.0;
t362 = t2.*t4.*t177.*v11.*v21.*7.0;
t363 = t2.*t7.*t177.*v11.*v21.*7.0;
t364 = t2.*t4.*t177.*v12.*v22.*7.0;
t365 = t4.*t7.*t177.*v12.*v22.*7.0;
t366 = t2.*t7.*t177.*v13.*v23.*7.0;
t367 = t4.*t7.*t177.*v13.*v23.*7.0;
t377 = t36+t178-2.0;
t378 = -r1.*t19.*(t15-t23);
t379 = -r2.*t19.*(t15-t23);
t380 = -r3.*t19.*(t15-t23);
t384 = r1.*r2.*r3.*t20.*(t15-t23).*-3.0;
t48 = -t41;
t49 = -t42;
t50 = -t43;
t54 = t41.*v11;
t58 = t42.*v12;
t62 = t43.*v13;
t69 = r2.*t55;
t73 = r3.*t56;
t74 = r3.*t59;
t203 = r3.*t171.*v11.*v21.*5.0;
t204 = r3.*t171.*v12.*v22.*5.0;
t205 = r3.*t171.*v13.*v23.*5.0;
t219 = t171.*v11.*v23.*-3.0;
t220 = t172.*v11.*v22.*-3.0;
t221 = t173.*v11.*v21.*-3.0;
t222 = t171.*v12.*v23.*-3.0;
t223 = t172.*v12.*v22.*-3.0;
t224 = t173.*v12.*v21.*-3.0;
t225 = t171.*v13.*v23.*-3.0;
t226 = t172.*v13.*v22.*-3.0;
t227 = t173.*v13.*v21.*-3.0;
t240 = t2.*t173.*v11.*v21;
t241 = t7.*t171.*v11.*v22;
t242 = t7.*t171.*v12.*v21;
t243 = t4.*t172.*v11.*v22;
t244 = t2.*t173.*v12.*v21;
t245 = t7.*t171.*v11.*v23;
t246 = t4.*t172.*v11.*v23;
t247 = t4.*t172.*v12.*v22;
t248 = t4.*t172.*v13.*v21;
t249 = t2.*t173.*v13.*v21;
t250 = t7.*t171.*v12.*v23;
t251 = t4.*t172.*v13.*v22;
t252 = t2.*t173.*v12.*v23;
t253 = t2.*t173.*v13.*v22;
t254 = t7.*t171.*v13.*v23;
t280 = -t206;
t281 = -t216;
t282 = -t217;
t283 = -t218;
t294 = r2.*t288;
t295 = r3.*t290;
t332 = r1.*t269.*v11.*7.0;
t333 = r2.*t268.*v12.*7.0;
t334 = r1.*t270.*v11.*7.0;
t335 = r3.*t268.*v13.*7.0;
t336 = r2.*t270.*v12.*7.0;
t337 = r3.*t269.*v13.*7.0;
t348 = r1.*t311.*v11.*7.0;
t349 = r2.*t309.*v12.*7.0;
t353 = r1.*t313.*v11.*7.0;
t355 = r3.*t308.*v13.*7.0;
t359 = r2.*t312.*v12.*7.0;
t360 = r3.*t310.*v13.*7.0;
t368 = -t347;
t369 = -t350;
t370 = -t351;
t371 = -t352;
t372 = -t354;
t373 = -t356;
t374 = -t357;
t375 = -t358;
t376 = -t361;
t385 = t27+t37+t177;
t386 = t18.*t377;
t299 = -t240;
t300 = -t243;
t301 = -t244;
t302 = -t245;
t303 = -t247;
t304 = -t249;
t305 = -t250;
t306 = -t251;
t307 = -t254;
et1 = t52+t53+t55-t56-t59+t61-t85+t86-t87+t88+t89+t92-t97+t98+t99+t100-t101-t102+t103-t105+t106-t108+t112+t117-t123-t124-t126+t127-t128+t129+t132+t133-t134-t136-t148-t153+t159+t160-t162-t163+t164+t166+t175+t176-t180+t181+t184-t186+t197+t199+t204-t205+t208+t210-t228+t229-t230+t232+t241+t246-t252-t253+t256+t257+t262+t263-t264-t266+t273+t276-t294+t295-t315+t316-t317+t318+t319+t322+t323-t324+t325-t327+t332-t333+t334-t335+t339+t342+t345+t348+t353-t359;
et2 = -t360+t362+t363-t364-t366+r3.*t42.*v13+t171.*v11.*v22.*6.0-t171.*v12.*v21.*9.0+t172.*v11.*v23.*6.0-t172.*v13.*v21.*9.0-t173.*v12.*v23.*3.0-t173.*v13.*v22.*3.0+r2.*t11.*t93.*v12+r3.*t11.*t93.*v13-r2.*t24.*v11.*v22.*4.0+r2.*t24.*v12.*v21.*6.0-r3.*t24.*v11.*v23.*4.0+r3.*t24.*v13.*v21.*6.0+r3.*t25.*v12.*v23.*2.0-t22.*t24.*v12.*v23.*2.0-t22.*t26.*v12.*v21.*2.0-r2.*t22.*t24.*v11.*v22.*2.0+r2.*t22.*t24.*v12.*v21.*3.0-r3.*t22.*t24.*v11.*v23.*2.0+r3.*t22.*t24.*v13.*v21.*3.0+r3.*t22.*t25.*v12.*v23+r3.*t22.*t25.*v13.*v22-r2.*t2.*t177.*v13.*v21.*8.0-r1.*t4.*t177.*v13.*v22.*8.0-t2.*t11.*t23.*v12.*v22.*3.0-t2.*t11.*t23.*v13.*v23.*3.0-r1.*r2.*r3.*t177.*v13.*v23.*8.0;
et3 = t72+t73+t74+t143+t145+t146+t156+t158+t161+t187+t225+t226+t227+t279+t290+t291+t292+t304+t306+t307+t373+t375+t376+r3.*t55-r3.*t175-t269.*v13.*8.0-t308.*v13.*8.0-t313.*v13.*8.0+r1.*t64.*v12.*1.5e+1-r2.*t65.*v11.*1.5e+1-r1.*t269.*v12.*7.0+r2.*t270.*v11.*7.0+r2.*t309.*v11.*7.0-r1.*t313.*v12.*7.0-t11.*t64.*v13.*5.0+t17.*v11.*v22.*2.0-t17.*v12.*v21-t24.*v13.*v21.*2.0+t25.*v12.*v23.*2.0-t25.*v13.*v22.*6.0+t26.*v12.*v22.*2.0-t26.*v13.*v23.*2.0+t171.*v11.*v21.*6.0-t171.*v12.*v22.*9.0-t172.*v12.*v23.*3.0+t173.*v11.*v23.*1.2e+1+r2.*t11.*t65.*v11+r1.*t11.*t94.*v12-r2.*t24.*v11.*v21.*4.0+r2.*t24.*v12.*v22.*6.0-r3.*t25.*v11.*v23.*8.0+r3.*t171.*v12.*v21.*5.0+t2.*t16.*v12.*v21.*2.0;
et4 = t4.*t16.*v11.*v22.*-1.0e+1+t4.*t16.*v12.*v21.*2.0+t6.*t15.*v13.*v22.*8.0-t7.*t16.*v11.*v22.*2.0+t8.*t15.*v11.*v22.*8.0-t8.*t23.*v11.*v22.*1.5e+1-t12.*t23.*v11.*v22.*2.0+t12.*t23.*v12.*v21-t22.*t25.*v12.*v23.*2.0+t22.*t25.*v13.*v22.*6.0-t22.*t26.*v12.*v22.*2.0+t24.*t36.*v13.*v21-t4.*t168.*v13.*v21.*5.0+t6.*t168.*v11.*v21+t4.*t170.*v12.*v22.*5.0-t4.*t170.*v13.*v23.*5.0+t6.*t170.*v11.*v23+t7.*t171.*v11.*v21-t4.*t172.*v12.*v23+t7.*t169.*v12.*v23.*5.0+t8.*t177.*v11.*v22.*7.0+r1.*t4.*t15.*v13.*v21.*8.0+r1.*t6.*t15.*v11.*v21.*8.0-r1.*t6.*t15.*v12.*v22.*8.0-r3.*t4.*t15.*v12.*v22.*8.0+r3.*t4.*t15.*v13.*v23.*8.0+r3.*t6.*t15.*v11.*v23.*8.0-r2.*t7.*t15.*v12.*v23.*8.0;
et5 = r2.*t9.*t15.*v11.*v23.*8.0-r1.*t6.*t23.*v11.*v21.*1.5e+1-r3.*t6.*t23.*v11.*v23.*1.5e+1+r1.*t12.*t23.*v13.*v21-r2.*t12.*t23.*v12.*v23+r2.*t12.*t31.*v13.*v22-r2.*t22.*t24.*v11.*v21.*2.0+r2.*t22.*t24.*v12.*v22.*3.0+r3.*t22.*t24.*v12.*v23-r3.*t22.*t25.*v11.*v23.*4.0+r1.*t6.*t177.*v11.*v21.*7.0+r3.*t4.*t177.*v12.*v22.*8.0+r3.*t6.*t177.*v11.*v23.*7.0+r2.*t7.*t177.*v12.*v23.*8.0-t2.*t4.*t15.*v12.*v21.*8.0+t4.*t7.*t15.*v11.*v22.*8.0+t2.*t4.*t23.*v12.*v21.*1.5e+1-t4.*t7.*t23.*v11.*v22.*1.5e+1-t2.*t11.*t23.*v12.*v21.*3.0+t4.*t11.*t23.*v11.*v22.*1.5e+1-t4.*t11.*t23.*v12.*v21.*3.0+t2.*t16.*t22.*v12.*v21-t4.*t16.*t22.*v11.*v22.*5.0+t4.*t16.*t22.*v12.*v21+t8.*t11.*t23.*v11.*v22-t7.*t16.*t22.*v11.*v22;
et6 = t6.*t16.*t22.*v13.*v22+t7.*t11.*t31.*v11.*v22+t4.*t22.*t24.*v13.*v21-t4.*t22.*t26.*v12.*v22+t4.*t22.*t26.*v13.*v23-t7.*t22.*t25.*v12.*v23-t2.*t4.*t177.*v12.*v21.*7.0+t4.*t7.*t177.*v11.*v22.*7.0-r1.*r2.*r3.*t15.*v12.*v21.*8.0+r1.*r2.*r3.*t177.*v12.*v21.*8.0+r1.*r2.*t7.*t15.*v11.*v21.*8.0-r1.*r3.*t4.*t15.*v12.*v23.*8.0-r1.*r2.*t7.*t23.*v11.*v21.*1.5e+1+r1.*r3.*t4.*t23.*v12.*v23.*1.5e+1-r2.*r3.*t22.*t24.*v12.*v21-t2.*t4.*t11.*t23.*v12.*v21+t4.*t7.*t11.*t23.*v11.*v22;
et7 = t58+t69+t70+t71+t140+t142+t144+t151+t154+t157-t183+t222+t223+t224+t287+t288+t289+t301+t303+t305+t370+t372+t374+r2.*t56+r2.*t176+t41.*v12+t270.*v12.*8.0+t309.*v12.*8.0+t311.*v12.*8.0-r3.*t64.*v11.*1.5e+1+r1.*t65.*v13.*1.5e+1+r3.*t269.*v11.*7.0-r1.*t270.*v13.*7.0+r3.*t308.*v11.*7.0-r1.*t311.*v13.*7.0+t11.*t65.*v12.*5.0+t17.*v11.*v23.*2.0-t17.*v13.*v21-t25.*v13.*v23.*2.0+t26.*v12.*v23.*6.0-t26.*v13.*v22.*2.0+t172.*v11.*v21.*6.0-t171.*v13.*v22.*3.0+t173.*v11.*v22.*1.2e+1-t172.*v13.*v23.*9.0+r3.*t11.*t64.*v11+r1.*t11.*t95.*v13-r3.*t24.*v11.*v21.*4.0-r3.*t25.*v11.*v22.*8.0+r3.*t24.*v13.*v23.*6.0-r3.*t171.*v13.*v21.*5.0+t2.*t16.*v13.*v21.*2.0-t4.*t16.*v11.*v23.*2.0;
et8 = t7.*t16.*v11.*v23.*-1.0e+1+t7.*t16.*v13.*v21.*2.0-t9.*t15.*v12.*v23.*8.0+t10.*t15.*v11.*v23.*8.0-t10.*t23.*v11.*v23.*1.5e+1-t12.*t23.*v11.*v23.*2.0+t12.*t23.*v13.*v21-t22.*t24.*v12.*v21.*2.0-t22.*t25.*v12.*v22.*2.0-t22.*t26.*v12.*v23.*6.0+t25.*t36.*v13.*v23+t26.*t36.*v13.*v22+t4.*t172.*v11.*v21+t7.*t168.*v12.*v21.*5.0-t4.*t170.*v13.*v22.*5.0+t9.*t168.*v11.*v21+t7.*t169.*v12.*v22.*5.0+t9.*t169.*v11.*v22-t7.*t169.*v13.*v23.*5.0-t7.*t171.*v13.*v22+t10.*t177.*v11.*v23.*7.0-r1.*t7.*t15.*v12.*v21.*8.0+r1.*t9.*t15.*v11.*v21.*8.0+r3.*t4.*t15.*v13.*v22.*8.0+r3.*t6.*t15.*v11.*v22.*8.0-r2.*t7.*t15.*v12.*v22.*8.0+r2.*t9.*t15.*v11.*v22.*8.0+r2.*t7.*t15.*v13.*v23.*8.0;
et9 = r1.*t9.*t15.*v13.*v23.*-8.0-r1.*t9.*t23.*v11.*v21.*1.5e+1-r2.*t9.*t23.*v11.*v22.*1.5e+1-r1.*t12.*t23.*v12.*v21-r3.*t12.*t23.*v12.*v23.*3.0+r3.*t12.*t23.*v13.*v22-r3.*t22.*t24.*v11.*v21.*2.0+r2.*t22.*t24.*v13.*v22-r3.*t22.*t25.*v11.*v22.*4.0+r3.*t22.*t24.*v13.*v23.*3.0+r1.*t9.*t177.*v11.*v21.*7.0-r3.*t4.*t177.*v13.*v22.*8.0+r2.*t9.*t177.*v11.*v22.*7.0-r2.*t7.*t177.*v13.*v23.*8.0-t2.*t7.*t15.*v13.*v21.*8.0+t4.*t7.*t15.*v11.*v23.*8.0+t2.*t7.*t23.*v13.*v21.*1.5e+1-t4.*t7.*t23.*v11.*v23.*1.5e+1-t2.*t11.*t23.*v13.*v21.*3.0+t2.*t16.*t22.*v13.*v21+t7.*t11.*t23.*v11.*v23.*1.5e+1-t7.*t11.*t23.*v13.*v21.*3.0-t4.*t16.*t22.*v11.*v23+t10.*t11.*t23.*v11.*v23-t7.*t16.*t22.*v11.*v23.*5.0;
et10 = t7.*t16.*t22.*v13.*v21+t4.*t11.*t31.*v11.*v23-t9.*t16.*t22.*v12.*v23-t7.*t22.*t24.*v12.*v21+t4.*t22.*t26.*v13.*v22-t7.*t22.*t25.*v12.*v22+t7.*t22.*t25.*v13.*v23-t2.*t7.*t177.*v13.*v21.*7.0+t4.*t7.*t177.*v11.*v23.*7.0+r1.*r2.*r3.*t15.*v13.*v21.*8.0-r1.*r2.*r3.*t177.*v13.*v21.*8.0+r1.*r3.*t4.*t15.*v11.*v21.*8.0-r1.*r2.*t7.*t15.*v13.*v22.*8.0-r1.*r3.*t4.*t23.*v11.*v21.*1.5e+1+r1.*r2.*t7.*t23.*v13.*v22.*1.5e+1+r2.*r3.*t22.*t24.*v13.*v21-t2.*t7.*t11.*t23.*v13.*v21+t4.*t7.*t11.*t23.*v11.*v23;
et11 = t62+t72+t73+t74+t143+t145+t146+t156+t158+t161-t187+t225+t226+t227+t290+t291+t292+t304+t306+t307+t373+t375+t376+r3.*t57+r3.*t174+r2.*t285+t42.*v13+t268.*v13.*8.0+t310.*v13.*8.0+t312.*v13.*8.0+r2.*t63.*v11.*1.5e+1-r1.*t65.*v12.*1.5e+1-r2.*t268.*v11.*7.0+r1.*t270.*v12.*7.0+r1.*t311.*v12.*7.0-r2.*t312.*v11.*7.0+t11.*t63.*v13.*5.0-t17.*v11.*v22+t17.*v12.*v21.*2.0-t24.*v11.*v23.*2.0+t24.*v13.*v21.*6.0-t26.*v11.*v21.*2.0-t171.*v11.*v21.*9.0+t171.*v12.*v22.*6.0+t172.*v12.*v23.*1.2e+1-t173.*v11.*v23.*3.0+r1.*t11.*t65.*v12+r2.*t11.*t93.*v11+r2.*t24.*v11.*v21.*6.0-r2.*t24.*v12.*v22.*4.0-r3.*t24.*v12.*v23.*8.0-r3.*t171.*v11.*v22.*5.0+t2.*t16.*v11.*v22.*2.0-t2.*t16.*v12.*v21.*1.0e+1;
et12 = t3.*t15.*v13.*v21.*-8.0+t4.*t16.*v11.*v22.*2.0+t5.*t15.*v12.*v21.*8.0-t7.*t16.*v12.*v21.*2.0-t5.*t23.*v12.*v21.*1.5e+1+t12.*t23.*v11.*v22-t12.*t23.*v12.*v21.*2.0-t22.*t24.*v13.*v21.*6.0-t22.*t25.*v13.*v22.*2.0-t22.*t26.*v13.*v23.*2.0+t24.*t36.*v11.*v23+t26.*t36.*v11.*v21-t2.*t170.*v11.*v21.*5.0+t2.*t169.*v13.*v22.*5.0+t3.*t169.*v12.*v22+t2.*t170.*v13.*v23.*5.0+t3.*t170.*v12.*v23-t2.*t173.*v11.*v23-t7.*t168.*v11.*v23.*5.0+t7.*t171.*v12.*v22+t5.*t177.*v12.*v21.*7.0-r2.*t3.*t15.*v11.*v21.*8.0+r3.*t2.*t15.*v11.*v21.*8.0-r2.*t2.*t15.*v13.*v22.*8.0+r2.*t3.*t15.*v12.*v22.*8.0-r3.*t2.*t15.*v13.*v23.*8.0+r3.*t3.*t15.*v12.*v23.*8.0+r1.*t7.*t15.*v11.*v23.*8.0+r1.*t9.*t15.*v12.*v23.*8.0;
et13 = r2.*t3.*t23.*v12.*v22.*-1.5e+1-r3.*t3.*t23.*v12.*v23.*1.5e+1+r1.*t12.*t23.*v11.*v23-r1.*t12.*t23.*v13.*v21.*3.0-r2.*t12.*t23.*v13.*v22+r2.*t22.*t24.*v11.*v21.*3.0-r2.*t22.*t24.*v12.*v22.*2.0-r3.*t22.*t24.*v12.*v23.*4.0+r3.*t22.*t25.*v11.*v23-r3.*t2.*t177.*v11.*v21.*8.0+r2.*t3.*t177.*v12.*v22.*7.0+r3.*t3.*t177.*v12.*v23.*7.0-r1.*t7.*t177.*v11.*v23.*8.0-t2.*t4.*t15.*v11.*v22.*8.0+t2.*t7.*t15.*v12.*v21.*8.0+t2.*t4.*t23.*v11.*v22.*1.5e+1-t2.*t7.*t23.*v12.*v21.*1.5e+1-t2.*t11.*t23.*v11.*v22.*3.0+t2.*t11.*t23.*v12.*v21.*1.5e+1-t4.*t11.*t23.*v11.*v22.*3.0+t5.*t11.*t23.*v12.*v21+t2.*t16.*t22.*v11.*v22-t2.*t16.*t22.*v12.*v21.*5.0-t3.*t16.*t22.*v13.*v21+t4.*t16.*t22.*v11.*v22;
et14 = -t7.*t16.*t22.*v12.*v21+t2.*t22.*t26.*v11.*v21+t7.*t11.*t31.*v12.*v21-t2.*t22.*t25.*v13.*v22-t2.*t22.*t26.*v13.*v23+t7.*t22.*t24.*v11.*v23-t2.*t4.*t177.*v11.*v22.*7.0+t2.*t7.*t177.*v12.*v21.*7.0+r1.*r2.*r3.*t15.*v11.*v22.*8.0-r1.*r2.*r3.*t177.*v11.*v22.*8.0-r2.*r3.*t2.*t15.*v11.*v23.*8.0+r1.*r2.*t7.*t15.*v12.*v22.*8.0+r2.*r3.*t2.*t23.*v11.*v23.*1.5e+1-r1.*r2.*t7.*t23.*v12.*v22.*1.5e+1-t2.*t4.*t11.*t23.*v11.*v22+t2.*t7.*t11.*t23.*v12.*v21;
et15 = t51+t53+t56-t57+t59-t60+t84+t85-t86-t90+t91+t92+t96-t98-t99-t100+t102-t103+t104+t107+t109-t110+t113-t117-t119+t123+t126-t127-t130+t131-t132+t134+t135-t137-t149+t153+t155-t159+t162-t164-t165+t167+t174+t176-t181+t182-t184+t185+t198+t201-t203+t205+t209+t213+t228+t230-t231-t233+t242-t246-t248+t252+t255+t257-t262+t264+t265-t267+t274+t277+t293-t295+t314+t315-t316-t320+t321+t322-t323-t325+t326+t328-t332+t333+t336-t337+t340+t341+t343+t349-t353;
et16 = -t355+t359-t362+t364+t365-t367+r3.*t41.*v13-t171.*v11.*v22.*9.0+t171.*v12.*v21.*6.0-t172.*v11.*v23.*3.0-t172.*v13.*v21.*3.0+t173.*v12.*v23.*6.0-t173.*v13.*v22.*9.0+r1.*t11.*t94.*v11+r3.*t11.*t94.*v13+r2.*t24.*v11.*v22.*6.0-r2.*t24.*v12.*v21.*4.0+r3.*t24.*v11.*v23.*2.0-r3.*t25.*v12.*v23.*4.0+r3.*t25.*v13.*v22.*6.0-t22.*t24.*v13.*v22.*2.0-t22.*t25.*v13.*v21.*2.0+r2.*t22.*t24.*v11.*v22.*3.0-r2.*t22.*t24.*v12.*v21.*2.0+r3.*t22.*t24.*v11.*v23+r3.*t22.*t24.*v13.*v21-r3.*t22.*t25.*v12.*v23.*2.0+r3.*t22.*t25.*v13.*v22.*3.0-r3.*t4.*t177.*v11.*v22.*8.0-r2.*t7.*t177.*v11.*v23.*8.0-t4.*t11.*t23.*v11.*v21.*3.0-t4.*t11.*t23.*v13.*v23.*3.0-r1.*r2.*r3.*t177.*v11.*v21.*8.0;
et17 = t66+t67+t68+t138+t139+t141+t147+t150+t152+t179+t219+t220+t221+t271+t284+t285+t286+t299+t300+t302+t368+t369+t371-r1.*t176-r2.*t291-t270.*v11.*8.0-t309.*v11.*8.0-t311.*v11.*8.0+r2.*t41.*v13-r3.*t63.*v12.*1.5e+1+r2.*t65.*v13.*1.5e+1+r3.*t268.*v12.*7.0-r2.*t270.*v13.*7.0-r2.*t309.*v13.*7.0+r3.*t310.*v12.*7.0-t11.*t65.*v11.*5.0+t17.*v12.*v23.*2.0-t17.*v13.*v22-t24.*v11.*v21.*2.0-t25.*v11.*v22.*2.0+t24.*v13.*v23.*2.0-t26.*v11.*v23.*6.0+t26.*v13.*v21.*2.0-t171.*v13.*v21.*3.0+t172.*v12.*v21.*1.2e+1+t173.*v12.*v22.*6.0-t173.*v13.*v23.*9.0+r3.*t11.*t63.*v12+r2.*t11.*t95.*v13-r3.*t24.*v12.*v21.*8.0-r3.*t25.*v12.*v22.*4.0+r3.*t25.*v13.*v23.*6.0+r3.*t171.*v13.*v22.*5.0;
et18 = t2.*t16.*v12.*v23.*-2.0+t4.*t16.*v13.*v22.*2.0-t7.*t16.*v12.*v23.*1.0e+1+t7.*t16.*v13.*v22.*2.0+t9.*t15.*v11.*v23.*8.0+t10.*t15.*v12.*v23.*8.0-t10.*t23.*v12.*v23.*1.5e+1-t12.*t23.*v12.*v23.*2.0+t12.*t23.*v13.*v22-t22.*t24.*v13.*v23.*2.0+t22.*t26.*v11.*v23.*6.0-t22.*t26.*v13.*v21.*2.0+t25.*t36.*v11.*v22+t2.*t170.*v13.*v21.*5.0-t7.*t168.*v11.*v21.*5.0+t2.*t173.*v12.*v22-t7.*t169.*v11.*v22.*5.0+t9.*t168.*v12.*v21+t7.*t168.*v13.*v23.*5.0-t7.*t171.*v13.*v21+t9.*t169.*v12.*v22+t10.*t177.*v12.*v23.*7.0-r3.*t2.*t15.*v13.*v21.*8.0+r3.*t3.*t15.*v12.*v21.*8.0+r1.*t7.*t15.*v11.*v21.*8.0+r2.*t7.*t15.*v11.*v22.*8.0+r1.*t9.*t15.*v12.*v21.*8.0-r1.*t7.*t15.*v13.*v23.*8.0+r2.*t9.*t15.*v12.*v22.*8.0;
et19 = r2.*t9.*t15.*v13.*v23.*-8.0-r1.*t9.*t23.*v12.*v21.*1.5e+1-r2.*t9.*t23.*v12.*v22.*1.5e+1+r2.*t12.*t23.*v11.*v22-r3.*t12.*t23.*v13.*v21+r3.*t12.*t31.*v11.*v23+r2.*t22.*t24.*v13.*v21-r3.*t22.*t24.*v12.*v21.*4.0-r3.*t22.*t25.*v12.*v22.*2.0+r3.*t22.*t25.*v13.*v23.*3.0+r3.*t2.*t177.*v13.*v21.*8.0+r1.*t9.*t177.*v12.*v21.*7.0+r1.*t7.*t177.*v13.*v23.*8.0+r2.*t9.*t177.*v12.*v22.*7.0+t2.*t7.*t15.*v12.*v23.*8.0-t4.*t7.*t15.*v13.*v22.*8.0-t2.*t7.*t23.*v12.*v23.*1.5e+1+t4.*t7.*t23.*v13.*v22.*1.5e+1-t4.*t11.*t23.*v13.*v22.*3.0-t2.*t16.*t22.*v12.*v23+t7.*t11.*t23.*v12.*v23.*1.5e+1-t7.*t11.*t23.*v13.*v22.*3.0+t4.*t16.*t22.*v13.*v22+t2.*t11.*t31.*v12.*v23+t10.*t11.*t23.*v12.*v23-t7.*t16.*t22.*v12.*v23.*5.0;
et20 = t7.*t16.*t22.*v13.*v22+t9.*t16.*t22.*v11.*v23-t2.*t22.*t26.*v13.*v21+t7.*t22.*t24.*v11.*v21+t7.*t22.*t25.*v11.*v22-t7.*t22.*t24.*v13.*v23+t2.*t7.*t177.*v12.*v23.*7.0-t4.*t7.*t177.*v13.*v22.*7.0-r1.*r2.*r3.*t15.*v13.*v22.*8.0+r1.*r2.*r3.*t177.*v13.*v22.*8.0+r2.*r3.*t2.*t15.*v12.*v22.*8.0-r1.*r2.*t7.*t15.*v13.*v21.*8.0-r2.*r3.*t2.*t23.*v12.*v22.*1.5e+1+r1.*r2.*t7.*t23.*v13.*v21.*1.5e+1+t2.*t7.*t11.*t23.*v12.*v23-t4.*t7.*t11.*t23.*v13.*v22;
et21 = t69+t70+t71+t140+t142+t144+t151+t154+t157+t183+t222+t223+t224+t275+t287+t288+t289+t301+t303+t305+t370+t372+t374-r2.*t174-r3.*t284-t268.*v12.*8.0-t310.*v12.*8.0-t312.*v12.*8.0+r3.*t42.*v11+r3.*t63.*v11.*1.5e+1-r1.*t64.*v13.*1.5e+1-r3.*t268.*v11.*7.0+r1.*t269.*v13.*7.0-r3.*t310.*v11.*7.0+r1.*t313.*v13.*7.0-t11.*t63.*v12.*5.0-t17.*v11.*v23+t17.*v13.*v21.*2.0+t24.*v11.*v22.*2.0-t24.*v12.*v21.*6.0+t25.*v11.*v21.*2.0-t25.*v12.*v22.*2.0-t26.*v12.*v23.*2.0-t172.*v11.*v21.*9.0+t171.*v13.*v22.*1.2e+1-t173.*v11.*v22.*3.0+t172.*v13.*v23.*6.0+r1.*t11.*t64.*v13+r3.*t11.*t93.*v11+r3.*t24.*v11.*v21.*6.0-r2.*t24.*v13.*v22.*8.0-r3.*t24.*v13.*v23.*4.0+r3.*t171.*v11.*v23.*5.0;
et22 = t3.*t15.*v12.*v21.*8.0+t2.*t16.*v11.*v23.*2.0-t2.*t16.*v13.*v21.*1.0e+1-t4.*t16.*v13.*v21.*2.0+t5.*t15.*v13.*v21.*8.0+t7.*t16.*v11.*v23.*2.0-t5.*t23.*v13.*v21.*1.5e+1+t12.*t23.*v11.*v23-t12.*t23.*v13.*v21.*2.0-t22.*t24.*v11.*v22.*2.0+t22.*t24.*v12.*v21.*6.0-t22.*t25.*v11.*v21.*2.0+t26.*t36.*v12.*v23+t2.*t169.*v11.*v21.*5.0-t2.*t169.*v12.*v22.*5.0+t4.*t168.*v11.*v22.*5.0-t2.*t170.*v12.*v23.*5.0+t3.*t169.*v13.*v22-t2.*t173.*v11.*v22+t3.*t170.*v13.*v23+t4.*t172.*v13.*v23+t5.*t177.*v13.*v21.*7.0-r2.*t2.*t15.*v11.*v21.*8.0-r1.*t4.*t15.*v11.*v22.*8.0+r2.*t2.*t15.*v12.*v22.*8.0-r3.*t3.*t15.*v11.*v21.*8.0+r2.*t3.*t15.*v13.*v22.*8.0+r3.*t2.*t15.*v12.*v23.*8.0+r1.*t6.*t15.*v13.*v22.*8.0;
et23 = r3.*t3.*t15.*v13.*v23.*8.0-r2.*t3.*t23.*v13.*v22.*1.5e+1-r3.*t3.*t23.*v13.*v23.*1.5e+1-r1.*t12.*t23.*v11.*v22+r3.*t12.*t23.*v12.*v23+r1.*t12.*t31.*v12.*v21+r3.*t22.*t24.*v11.*v21.*3.0-r2.*t22.*t24.*v13.*v22.*4.0+r3.*t22.*t25.*v11.*v22-r3.*t22.*t24.*v13.*v23.*2.0+r2.*t2.*t177.*v11.*v21.*8.0+r1.*t4.*t177.*v11.*v22.*8.0+r2.*t3.*t177.*v13.*v22.*7.0+r3.*t3.*t177.*v13.*v23.*7.0+t2.*t4.*t15.*v13.*v21.*8.0-t2.*t7.*t15.*v11.*v23.*8.0-t2.*t4.*t23.*v13.*v21.*1.5e+1+t2.*t7.*t23.*v11.*v23.*1.5e+1-t2.*t11.*t23.*v11.*v23.*3.0+t2.*t11.*t23.*v13.*v21.*1.5e+1+t5.*t11.*t23.*v13.*v21+t2.*t16.*t22.*v11.*v23-t2.*t16.*t22.*v13.*v21.*5.0+t3.*t16.*t22.*v12.*v21-t7.*t11.*t23.*v11.*v23.*3.0-t4.*t16.*t22.*v13.*v21;
et24 = t7.*t16.*t22.*v11.*v23+t4.*t11.*t31.*v13.*v21-t2.*t22.*t25.*v11.*v21+t2.*t22.*t25.*v12.*v22-t4.*t22.*t24.*v11.*v22+t2.*t22.*t26.*v12.*v23+t2.*t4.*t177.*v13.*v21.*7.0-t2.*t7.*t177.*v11.*v23.*7.0-r1.*r2.*r3.*t15.*v11.*v23.*8.0+r1.*r2.*r3.*t177.*v11.*v23.*8.0-r2.*r3.*t2.*t15.*v11.*v22.*8.0+r1.*r3.*t4.*t15.*v13.*v23.*8.0+r2.*r3.*t2.*t23.*v11.*v22.*1.5e+1-r1.*r3.*t4.*t23.*v13.*v23.*1.5e+1+t2.*t4.*t11.*t23.*v13.*v21-t2.*t7.*t11.*t23.*v11.*v23;
et25 = t54+t66+t67+t68+t138+t139+t141+t147+t150+t152-t179+t219+t220+t221+t284+t285+t286+t299+t300+t302+t368+t369+t371+r1.*t175+r3.*t287+t43.*v11+t269.*v11.*8.0+t308.*v11.*8.0+t313.*v11.*8.0+r3.*t41.*v12-r2.*t63.*v13.*1.5e+1+r3.*t64.*v12.*1.5e+1+r2.*t268.*v13.*7.0-r3.*t269.*v12.*7.0-r3.*t308.*v12.*7.0+r2.*t312.*v13.*7.0+t11.*t64.*v11.*5.0-t17.*v12.*v23+t17.*v13.*v22.*2.0-t24.*v12.*v22.*2.0+t25.*v11.*v22.*6.0-t25.*v12.*v21.*2.0+t171.*v13.*v21.*1.2e+1-t172.*v12.*v21.*3.0-t173.*v12.*v22.*9.0+t173.*v13.*v23.*6.0+r2.*t11.*t63.*v13+r3.*t11.*t94.*v12-r2.*t24.*v13.*v21.*8.0+r3.*t25.*v12.*v22.*6.0-r3.*t25.*v13.*v23.*4.0-r3.*t171.*v12.*v23.*5.0-t2.*t16.*v13.*v22.*2.0-t6.*t15.*v11.*v22.*8.0;
et26 = t4.*t16.*v12.*v23.*2.0-t4.*t16.*v13.*v22.*1.0e+1+t7.*t16.*v12.*v23.*2.0+t8.*t15.*v13.*v22.*8.0-t8.*t23.*v13.*v22.*1.5e+1+t12.*t23.*v12.*v23-t12.*t23.*v13.*v22.*2.0-t22.*t24.*v11.*v21.*2.0-t22.*t25.*v11.*v22.*6.0-t22.*t26.*v11.*v23.*2.0+t24.*t36.*v12.*v22+t25.*t36.*v12.*v21-t2.*t169.*v12.*v21.*5.0+t4.*t168.*v11.*v21.*5.0-t4.*t168.*v12.*v22.*5.0+t4.*t170.*v11.*v23.*5.0+t6.*t168.*v13.*v21-t4.*t172.*v12.*v21+t2.*t173.*v13.*v23+t6.*t170.*v13.*v23+t8.*t177.*v13.*v22.*7.0-r1.*t4.*t15.*v11.*v21.*8.0+r2.*t2.*t15.*v12.*v21.*8.0+r1.*t4.*t15.*v12.*v22.*8.0+r2.*t3.*t15.*v13.*v21.*8.0+r1.*t6.*t15.*v13.*v21.*8.0-r3.*t4.*t15.*v11.*v23.*8.0-r3.*t6.*t15.*v12.*v22.*8.0+r3.*t6.*t15.*v13.*v23.*8.0;
et27 = r1.*t6.*t23.*v13.*v21.*-1.5e+1-r3.*t6.*t23.*v13.*v23.*1.5e+1-r2.*t12.*t23.*v11.*v22.*3.0+r2.*t12.*t23.*v12.*v21-r3.*t12.*t23.*v11.*v23-r2.*t22.*t24.*v13.*v21.*4.0+r3.*t22.*t24.*v12.*v21+r3.*t22.*t25.*v12.*v22.*3.0-r3.*t22.*t25.*v13.*v23.*2.0-r2.*t2.*t177.*v12.*v21.*8.0-r1.*t4.*t177.*v12.*v22.*8.0+r1.*t6.*t177.*v13.*v21.*7.0+r3.*t6.*t177.*v13.*v23.*7.0+t2.*t4.*t15.*v13.*v22.*8.0-t4.*t7.*t15.*v12.*v23.*8.0-t2.*t4.*t23.*v13.*v22.*1.5e+1+t4.*t7.*t23.*v12.*v23.*1.5e+1-t4.*t11.*t23.*v12.*v23.*3.0+t4.*t11.*t23.*v13.*v22.*1.5e+1-t2.*t16.*t22.*v13.*v22-t7.*t11.*t23.*v12.*v23.*3.0+t4.*t16.*t22.*v12.*v23-t4.*t16.*t22.*v13.*v22.*5.0-t6.*t16.*t22.*v11.*v22+t8.*t11.*t23.*v13.*v22;
et28 = t2.*t11.*t31.*v13.*v22+t7.*t16.*t22.*v12.*v23+t2.*t22.*t25.*v12.*v21-t4.*t22.*t24.*v11.*v21+t4.*t22.*t24.*v12.*v22-t4.*t22.*t26.*v11.*v23+t2.*t4.*t177.*v13.*v22.*7.0-t4.*t7.*t177.*v12.*v23.*7.0+r1.*r2.*r3.*t15.*v12.*v23.*8.0-r1.*r2.*r3.*t177.*v12.*v23.*8.0-r1.*r3.*t4.*t15.*v12.*v21.*8.0+r2.*r3.*t2.*t15.*v13.*v23.*8.0+r1.*r3.*t4.*t23.*v12.*v21.*1.5e+1-r2.*r3.*t2.*t23.*v13.*v23.*1.5e+1+t2.*t4.*t11.*t23.*v13.*v22-t4.*t7.*t11.*t23.*v12.*v23;
et29 = t51+t52-t55+t57+t60-t61+t84+t87+t88-t89+t90-t91-t96+t97+t101-t104+t105-t106-t107+t108-t109+t110-t112-t113+t119+t124+t128-t129+t130-t131-t133-t135+t136+t137+t148+t149-t155-t160+t163+t165-t166-t167+t174+t175+t180-t182-t185+t186+t200+t202+t203-t204+t212+t214-t229+t231-t232+t233-t241-t242+t248+t253+t255+t256-t263-t265+t266+t267+t272+t278-t293+t294+t314+t317+t318-t319+t320-t321+t324-t326+t327-t328-t334+t335-t336+t337+t338+t344+t346-t348;
et30 = -t349+t355+t360-t363-t365+t366+t367+r2.*t41.*v12-t171.*v11.*v22.*3.0-t171.*v12.*v21.*3.0-t172.*v11.*v23.*9.0+t172.*v13.*v21.*6.0-t173.*v12.*v23.*9.0+t173.*v13.*v22.*6.0+r1.*t11.*t95.*v11+r2.*t11.*t95.*v12+r2.*t24.*v11.*v22.*2.0+r3.*t24.*v11.*v23.*6.0-r3.*t24.*v13.*v21.*4.0+r3.*t25.*v12.*v23.*6.0-r3.*t25.*v13.*v22.*4.0-t22.*t25.*v11.*v23.*2.0-t22.*t26.*v11.*v22.*2.0+r2.*t22.*t24.*v11.*v22+r2.*t22.*t24.*v12.*v21+r3.*t22.*t24.*v11.*v23.*3.0-r3.*t22.*t24.*v13.*v21.*2.0+r3.*t22.*t25.*v12.*v23.*3.0-r3.*t22.*t25.*v13.*v22.*2.0-r3.*t2.*t177.*v12.*v21.*8.0-r1.*t7.*t177.*v12.*v23.*8.0-t7.*t11.*t23.*v11.*v21.*3.0-t7.*t11.*t23.*v12.*v22.*3.0-r1.*r2.*r3.*t177.*v12.*v22.*8.0;
dAv1drTv2_dr = reshape([-t21.*(et1+et2),-t21.*(et11+et12+et13+et14),-t21.*(et21+et22+et23+et24),-t21.*(et3+et4+et5+et6),-t21.*(et15+et16),-t21.*(et25+et26+et27+et28),-t21.*(et7+et8+et9+et10),-t21.*(et17+et18+et19+et20),-t21.*(et29+et30)],[3,3]);
if nargout > 1
    t381 = t28+t281;
    t382 = t29+t282;
    t383 = t30+t283;
    t387 = -t386;
    t391 = r1.*t13.*t385;
    t392 = r2.*t13.*t385;
    t393 = r3.*t13.*t385;
    t388 = r2.*r3.*t19.*t381;
    t389 = r1.*r3.*t19.*t382;
    t390 = r1.*r2.*t19.*t383;
    mt1 = [t20.*(t49+t50+t76+t77+t78+t79+t81+t83+t94+t95+t269+t270+t308+t309+t311+t313-r1.*t4.*t23.*v21.*3.0+r2.*t2.*t23.*v22.*2.0+r3.*t2.*t23.*v23.*2.0-r1.*t7.*t23.*v21.*3.0-r3.*t4.*t23.*v23-r2.*t7.*t23.*v22),-v21.*(t238+t297+t379-r1.*r2.*t19.*t381+r2.*t2.*t20.*(t15-t23).*3.0)-v23.*(t261+t280+t331-t390+r1.*r2.*r3.*t20.*(t15-t23).*3.0)-t20.*v22.*(-t24+t47+t168+t173+t236-r1.*t4.*t23.*3.0+r1.*t4.*t27+r1.*t4.*t177)];
    mt2 = [v22.*(t260+t280+t330+t384+t389)+v21.*(t237+t296+r3.*t19.*(t15-t23)+r1.*r3.*t19.*t381-r3.*t2.*t20.*(t15-t23).*3.0)+t20.*v23.*(t24+t47-t168+t173+t236-r1.*t7.*t15.*2.0+r1.*t7.*t31+r1.*t7.*t22.*t35),v23.*(t261+t280+t331+t384+t390)+v22.*(t239+t298+r1.*t19.*(t15-t23)+r1.*r2.*t19.*t382-r1.*t4.*t20.*(t15-t23).*3.0)+t20.*v21.*(t25+t46-t169+t172+t235-r2.*t2.*t15.*2.0+r2.*t2.*t31+r2.*t2.*t22.*t35),t20.*(t48+t50+t75+t77+t79+t80+t81+t82+t93+t95+t268+t270+t309+t310+t311+t312+r1.*t4.*t23.*v21.*2.0-r2.*t2.*t23.*v22.*3.0-r3.*t2.*t23.*v23-r1.*t7.*t23.*v21+r3.*t4.*t23.*v23.*2.0-r2.*t7.*t23.*v22.*3.0)];
    mt3 = [-v22.*(t237+t296+t380-r2.*r3.*t19.*t382+r3.*t4.*t20.*(t15-t23).*3.0)-v21.*(t259+t280+t329-t388+r1.*r2.*r3.*t20.*(t15-t23).*3.0)-t20.*v23.*(-t25+t46+t169+t172+t235-r2.*t7.*t23.*3.0+r2.*t7.*t27+r2.*t7.*t177),-v23.*(t239+t298+t378-r1.*r3.*t19.*t383+r1.*t7.*t20.*(t15-t23).*3.0)-v22.*(t260+t280+t330-t389+r1.*r2.*r3.*t20.*(t15-t23).*3.0)-t20.*v21.*(-t26+t45+t170+t171+t234-r3.*t2.*t23.*3.0+r3.*t2.*t27+r3.*t2.*t177)];
    mt4 = [v21.*(t259+t280+t329+t384+t388)+v23.*(t238+t297+r2.*t19.*(t15-t23)+r2.*r3.*t19.*t383-r2.*t7.*t20.*(t15-t23).*3.0)+t20.*v22.*(t26+t45-t170+t171+t234-r3.*t4.*t15.*2.0+r3.*t4.*t31+r3.*t4.*t22.*t35),t20.*(t48+t49+t75+t76+t78+t80+t82+t83+t93+t94+t268+t269+t308+t310+t312+t313-r1.*t4.*t23.*v21-r2.*t2.*t23.*v22-r3.*t2.*t23.*v23.*3.0+r1.*t7.*t23.*v21.*2.0-r3.*t4.*t23.*v23.*3.0+r2.*t7.*t23.*v22.*2.0)];
    dAv1drTv2_dv1 = reshape([mt1,mt2,mt3,mt4],3,3);
end
if nargout > 2
    C = reshape([t387,-t393,t392,t393,t387,-t391,-t392,t391,t387],[3,3]);
end
