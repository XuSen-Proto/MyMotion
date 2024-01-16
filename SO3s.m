classdef SO3s < handle
    properties (Constant)
        tol = 1e-12;
    end
    properties
        rn = zeros(3,1);
        th = 0;
        th_1 = 0;
        Sth = 1;
        Cth = 0;
        R = eye(3);
        c0 = 0;
        c1 = 0;
        c2 = 0;
        c3 = 0;
        A = eye(3);
    end
    methods
        function obj = SO3s(Arg1,th)
            if size(Arg1,2)>1
                obj.LogUpd(Arg1);
            else
                obj.ExpUpd(Arg1,th);
            end
        end

        function dwb = dwbdso3E(obj,dr)
            obj.c1 = obj.c0*obj.th_1;
            obj.c2 = (obj.th-obj.Sth)*obj.th_1;
            obj.A = so3upoly([1 -obj.c1 obj.c2],obj.rn); 
            dwb = obj.A*dr;
        end

        function C = d_dwbdso3E_Pt(obj,t)
            t = t*obj.th_1;
            obj.c3 = 2*obj.c1-obj.Sth;
            rht = Ahat3B(obj.rn,t);
            Cthat = so3upoly([obj.c1 obj.c2-obj.c0 obj.c3],obj.rn);
            Cthat = Ahat3B(t,Cthat');
            CtrT = so3upoly([-obj.c3 obj.Sth*obj.c1-(2*obj.Cth+1)*obj.c2...
                obj.Sth*(2+obj.c2)-4*obj.c1],obj.rn);
            C = CtrT*rht*obj.rn'+Cthat'+Ahat3B(obj.c2*rht,obj.R)';
        end

        function dr = dso3dwbE(obj,dwb)
            obj.c1 = -obj.c0*obj.th_1;
            obj.c2 = obj.th/2;
            obj.c3 = 1+obj.Sth/(2*obj.c1);
            obj.A = so3upoly([1 obj.c2 obj.c3],obj.rn);
            dr = obj.A*dwb;
        end

        function C = d_dso3dwbE_Pt(obj,t)
            if obj.th~=0
                t = t*obj.th_1;
                rht = Ahat3B(obj.rn,t);
                obj.c3 = 1+obj.Sth/(2*obj.c1);
                Cthat = so3upoly([obj.c1 obj.c3 obj.c1+obj.c2],obj.rn);
                CtrT = so3upoly(-[obj.Sth+2*obj.c1 3*obj.c3-0.5*obj.c0...
                    2*obj.c1+obj.Sth/2+obj.c2],obj.rn);
                
                Cthat = Ahat3B(t,Cthat');
                C = (CtrT*rht*obj.rn'+Cthat'+...
                    Ahat3B((obj.Sth/obj.th-1)*rht,obj.A)')*obj.A;
            else
                C = eye(3);
            end
        end

        function dws = dwsdrE(obj,dr)
            obj.c1 = obj.c0*obj.th_1;
            obj.c2 = (obj.th-obj.Sth)*obj.th_1;
            obj.A = so3upoly([1 obj.c1 obj.c2],obj.rn); 
            dws = obj.A*dr;
        end

        function C = d_dwsdrE_Ptn(obj,t)
            rht = Ahat3B(obj.rn,t);
            Cthat = hat3(-obj.c2*obj.rn);
            Cthat(1:4:9) = Cthat(1:4:9)-obj.c1;
            Cthat = Ahat3B(obj.rn,Cthat);
            Cthat = Ahat3B(t,Cthat');
            CtrT = hat3((obj.c0-3*obj.c2)*obj.rn);
            CtrT(1:4:9) = CtrT(1:4:9)+(obj.Sth-2*obj.c1);
            CtrT = Ahat3B(obj.rn,CtrT);
            C = CtrT*rht*obj.rn'+Cthat'-hat3(obj.c2*rht)+Ahat3B(obj.p,obj.A);
        end 

        function dr = drdwsE(obj,dwb)
            obj.c1 = obj.c0*obj.th_1;
            obj.A = so3upoly([1 -obj.th/2 1-obj.Sth/(2*obj.c1)],obj.rn);
            dr = obj.A*dwb;
        end

        function C = d_drdwsE_Ptn(obj,t)
            obj.c2 = (obj.th-obj.Sth)*obj.th_1;
            rht = Ahat3B(obj.rn,t);
            Cthat = hat3(-obj.c2*obj.rn);
            Cthat(1:4:9) = Cthat(1:4:9)-obj.c1;
            Cthat = Ahat3B(obj.rn,Cthat);
            Cthat = Ahat3B(t,Cthat');
            CtrT = hat3((obj.c0-3*obj.c2)*obj.rn);
            CtrT(1:4:9) = CtrT(1:4:9)+(obj.Sth-2*obj.c1);
            CtrT = Ahat3B(obj.rn,CtrT);
            C = (CtrT*rht*obj.rn'+Cthat'-hat3(obj.c2*rht))*obj.A+hat3(obj.p);
            C = -obj.A*C;
        end 
        function LogUpd(obj,Rin)
            Cth = (Rin(1,1)+Rin(2,2)+Rin(3,3)-1)/2;
            Cth = min(max(Cth,-1),1);
            if (1-obj.Cth)<obj.tol
                obj.rn = zeros(3,1);
                obj.th = 0;
                obj.Sth = 0;
                obj.Cth = 1;
                obj.R = eye(3);
            elseif (obj.Cth+1)<obj.tol
                obj.th = pi;
                obj.Sth = 0;
                obj.Cth = -1;
                obj.rn = [Rin(1:2,3); Rin(3,3)+1];
                obj.rn = (1/sqrt(2*obj.rn(3)))*obj.rn;
                obj.R = Rin;
            else
                obj.R = Rin;
                obj.Cth = Cth;
                obj.th = acos(obj.Cth); 
                obj.Sth = sin(obj.th);
                obj.rn = (1/(2*obj.Sth))*[Rin(3,2)-Rin(2,3); ...
                    Rin(1,3)-Rin(3,1); Rin(2,1)-Rin(1,2)];
            end
            obj.c0 = 1-obj.Cth;
            obj.th_1 = 1/obj.th;
        end

        function ExpUpd(obj,omg,th)
            if isempty(th)
                th = norm(omg);
                if abs(th)<obj.tol
                    obj.th = 0;
                    obj.rn = zeros(3,1);
                    obj.th_1 = inf;
                    obj.Sth = 0;
                    obj.Cth = 1;
                    obj.c0 = 0;
                    obj.R = eye(3);
                    return;
                else
                    obj.th = th;
                    obj.th_1 = 1/th;
                    obj.rn = omg*obj.th_1;
                end
            else
                obj.th = th;
                obj.rn = omg;
            end
            obj.Sth = sin(th);
            obj.Cth = cos(th);
            obj.c0 = 1-obj.Cth;
            obj.R = so3upoly([1 obj.Sth obj.c0],obj.rn);
            obj.th_1 = 1/obj.th;
        end
    end
end