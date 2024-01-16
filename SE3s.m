classdef SE3s < handle
    properties (Constant)
        tol = 1e-12;
    end
    properties
        rn = zeros(3,1);
        vn = zeros(3,1);
        p = zeros(3,1);
        RTp = zeros(3,1);
        th = 0;
        th_1 = 0;
        Sth = 1;
        Cth = 0;
        R = eye(3);
        c0 = 0;
        A = eye(3); Avn = zeros(3); C = zeros(3);% analytical jacobian of SO3
        phR = zeros(3);
        dArhv = zeros(3,1);
        rhv = zeros(3,1);
        C1 = zeros(3);
        C2 = zeros(3);
        a = [0 0 0];
        da = [0 0];
        dda = [0 0];
    end
    methods
        function obj = SE3s(Arg1,th)
            if size(Arg1,2)>1
                obj.LogUpd(Arg1);
            else
                obj.ExpUpd(Arg1(1:3),Arg1(4:6),th);
            end
        end

        function GetRTp(obj)
            obj.RTp = (obj.R'*obj.p);
        end

        function GetphR(obj)
            obj.phR = hat3(obj.p)*obj.R;
        end

        function [c1,c2,c3]=dwdrCoefs(obj,Rtag)
            % Rtag: 0-ws, 1-wb
            sgn = (-1)^Rtag;
            ct1 = obj.c0/obj.th;
            ct2 = 1-obj.Sth/obj.th;
            c1 = [1 sgn*ct1 ct2];
            c2 = [sgn*(obj.Sth-2*ct1) obj.c0-3*ct2];
            c3 = [sgn*(obj.th*obj.Cth+8*ct1-5*obj.Sth) ...
                15*ct2-7*obj.c0+obj.th*obj.Sth];
        end

        function [c1,c2,c3]=drdwCoefs(obj,Rtag)
            ct1 = obj.th/2;
            ct2 = 1-ct1*obj.Sth/obj.c0;
            ct3 = ct1*(obj.th+obj.Sth)/obj.c0;
            c1 = [1 (-(-1)^Rtag)*ct1 ct2];
            c2 = [0 ct3-2];
            c3 = [0 8-3*ct3-obj.th^3*obj.Sth/(2*obj.c0^2)];
        end

        function getA(obj,c)
            obj.A = so3upoly(c,obj.rn);
        end

        function dA = dAvn_dr(obj,vin,c,cd)
            obj.rhv = Ahat3B(obj.rn,vin);
            dA = hat3(cd(2)*obj.rn);
            dA(1:4:9) = dA(1:4:9)+cd(1);
            obj.C1 = dA;
            obj.dArhv = dA*obj.rhv;
            dA = obj.dArhv*obj.rn';
            Cthat = hat3(-c(2)*obj.rn);
            Cthat(1:4:9) = Cthat(1:4:9)-c(1);
            Cthat = Cthat*hat3(vin);
            dA = dA + Cthat - c(2)*hat3(obj.rhv);
        end

        function [dAdr,dAdv] = ddAvndrTmn(obj,vin,mn,c,cd,cdd)
            rhm = Ahat3B(obj.rn,mn);
            dAdv = hat3(-c(2)*obj.rn);
            dAdv(1:4:9) = dAdv(1:4:9)-c(1);
            dAdv = Ahat3B(mn,dAdv);
            dAdv = dAdv+hat3(c(2)*rhm)-obj.rn*(rhm'*obj.C1);
            dAdr = hat3(-cd(2)*obj.rn);
            dAdr(1:4:9) = dAdr(1:4:9)+cd(1);
            dAdr = Ahat3B(vin,dAdr)+hat3(cd(2)*obj.rhv);
            dAdr = dAdr*mn*obj.rn';
            dAdr(1:4:9) = dAdr(1:4:9)+obj.dArhv'*mn;
            dAdr = dAdr+obj.rn*(vin'*obj.dAvn_dr(mn,[-cd(1) cd(2)],[-cdd(1) cdd(2)]));
            mn = c(2)*mn;
            dAdr = dAdr+hat3(mn)*hat3(vin)+hat3(vin)*hat3(mn);
        end

        function [dAdr,dAdv] = ddAvndrmn(obj,vin,mn,c,cd,cdd)
            vhm = Ahat3B(vin,mn);
            C1t = obj.C1*hat3(obj.rn);
            rTm = obj.rn'*mn;
            dAdv = hat3(c(2)*obj.rn);
            dAdv(1:4:9) = dAdv(1:4:9)+c(1);
            dAdv = Ahat3B(mn,dAdv);
            dAdr = hat3(-cd(2)*obj.rn);
            dAdr(1:4:9) = dAdr(1:4:9)-cd(1);
            dAdr = (dAdr*vhm-Ahat3B(obj.rhv,cd(2)*mn))*obj.rn'+obj.dArhv*mn';
            mn = c(2)*mn;
            dAdv = dAdv+hat3(obj.rn)*hat3(mn)+rTm*C1t;
            dAdr = dAdr+rTm*obj.dAvn_dr(vin,cd,cdd)+hat3(c(2)*vhm)-...
                Ahat3B(mn,hat3(vin));
        end

        function [A,C] = dVdr(obj,Rtag)
            [obj.a,obj.da,obj.dda]= obj.dwdrCoefs(Rtag);
            obj.getA(obj.a);
%             obj.GetRTp();
            A = obj.A;
            obj.Avn = obj.dAvn_dr(obj.vn,obj.a(2:3),obj.da);
            C = obj.Avn-Ahat3B(obj.RTp,A);
        end

        function [A,C] = drdV(obj,Rtag)
            [c1,c2,~]= obj.drdwCoefs(Rtag);
            [b1,b2,~]= obj.dwdrCoefs(Rtag);
            obj.getA(c1);
%             obj.GetRTp();
            A = obj.A;
            C = A*(hat3(obj.RTp)-obj.dAvn_dr(obj.vn,b1,b2)*A);
        end

        function LogUpd(obj,T)
            obj.R = T(1:3,1:3);
            obj.p = T(1:3,4);
            CthTmp = (obj.R(1,1)+obj.R(2,2)+obj.R(3,3)-1)/2;
            CthTmp = min(max(CthTmp,-1),1);
            c0Tmp = 1-CthTmp;
            cp = (CthTmp+1)<obj.tol;
            if (c0Tmp>=obj.tol)&&(~cp)
                obj.th = acos(CthTmp);
                obj.Cth = CthTmp;
                obj.c0 = c0Tmp;
                obj.Sth = sin(obj.th);
                obj.rn = (1/(2*obj.Sth))*[obj.R(3,2)-obj.R(2,3); ...
                    obj.R(1,3)-obj.R(3,1); obj.R(2,1)-obj.R(1,2)];
                obj.th_1 = 1/obj.th;
                obj.vn = so3upoly([obj.th_1 -1/2 obj.th_1-obj.Sth/(2*obj.c0)],...
                    obj.rn)*obj.p;
            elseif cp
                obj.th = pi;
                v3 = obj.R(3,3)+1;
                obj.rn = [obj.R(1:2,3); v3];
                obj.rn = (1/sqrt(2*v3))*obj.rn;
                obj.vn = so3upoly([obj.th_1 -1/2 obj.th_1],obj.rn)*obj.p;
                return;
            else
                obj.rn = zeros(3,1);
                obj.vn = obj.p;
                obj.th = 0;
                obj.Sth = 0;
                obj.Cth = 1;
                obj.c0 = 0;
                return;
            end
        end

        function ExpUpd(obj,Omg,v,th)
            if isempty(th)
                th = norm(Omg);
                if th<obj.tol
                    obj.th = 0;
                    obj.rn = zeros(3,1);
                    obj.th_1 = inf;
                    obj.Sth = 0;
                    obj.Cth = 1;
                    obj.c0 = 0;
                    obj.R = eye(3);
                    obj.vn = v;
                    obj.p = v;
                else
                    obj.th = th;
                    obj.Sth = sin(th);
                    obj.Cth = cos(th);
                    obj.rn = Omg/th;
                    obj.vn = v/th;
                    obj.c0 = 1-obj.Cth;
                    obj.R = so3upoly([1 obj.Sth obj.c0],obj.rn);
                    obj.p = so3upoly([th obj.c0 th-obj.Sth],obj.rn)*obj.vn;
                end
            else
                obj.th = th;
                obj.rn = Omg;
                obj.vn = v;
                obj.Sth = sin(th);
                obj.Cth = cos(th);
                obj.c0 = 1-obj.Cth;
                obj.R = so3upoly([1 obj.Sth obj.c0],obj.rn);
                obj.p = so3upoly([th obj.c0 th-obj.Sth],obj.rn)*obj.vn;
            end
        end

        function S = Ad(obj,S)
            S = [obj.R*S(1:3,:); obj.R*S(4:6,:)+obj.phR*S(1:3,:)];
        end

        function S = Adi(obj,S)
            S = [obj.R'*S(1:3,:); obj.R'*S(4:6,:)+obj.phR'*S(1:3,:)];
        end

        function S = AdT(obj,S)
            S = [obj.R'*S(1:3,:)+obj.phR'*S(4:6,:); obj.R'*S(4:6,:)];
        end

        function S = AdTi(obj,S)
            S = [obj.R*S(1:3,:)+obj.phR*S(4:6,:); obj.R*S(4:6,:)];
        end

        function S = AdIdx(obj,i)
            if i<4
                S = [obj.R(:,i); obj.phR(:,i)];
            else
                S = [zeros(3,1); obj.R(:,i-3)];
            end
        end

        function S = AdiIdx(obj,i)
            if i<4
                S = [obj.R(i,:)'; obj.phR(i,:)'];
            else
                S = [zeros(3,1); obj.R(i-3,:)'];
            end
        end

        function S = AdTIdx(obj,i)
            if i<4
                S = [obj.R(i,:)'; zeros(3,1)];
            else
                i = i-3;
                S = [obj.phR(i,:)'; obj.R(i,:)'];
            end
        end

        function S = AdTiIdx(obj,i)
            if i<4
                S = [obj.R(:,i); zeros(3,1)];
            else
                i = i-3;
                S = [obj.phR(:,i); obj.R(:,i)];
            end
        end

    end
end