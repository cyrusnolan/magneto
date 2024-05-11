classdef state_vec_N
    %STATE_VEC_N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NRao0
        NRbo0
        NRfo0
        NRgo0
        NRdo0
        NReo0
        NVao0
        NVbo0
        NVfo0
        NVgo0
        NVdo0
        NVeo0
    end
    
    methods
        function obj = state_vec_N(p, orbit)
            %STATE_VEC_N Construct an instance of this class
            %   Detailed explanation goes here
            
            rp = sqrt((p.l0-p.p0)^2 - p.d0^2);
            % form state vec in body coordinates
            b1 = [1;0;0]; % basis vectors
            b2 = [0;1;0];
            b3 = [0;0;1];
            BWc0 = p.w*b3;
            BRac0 = rp*b1;
            BRbc0 = -rp*b1;
            BRfc0 = rp*b2;
            BRgc0 = -rp*b2;
            BRdc0 = p.d0*b3;
            BRec0 = -p.d0*b3;
            BVac0 = cross(BWc0, BRac0);
            BVbc0 = cross(BWc0, BRbc0);
            BVfc0 = cross(BWc0, BRfc0);
            BVgc0 = cross(BWc0, BRgc0);
            BVdc0 = zeros(3, 1);
            BVec0 = zeros(3, 1);

            % convert state vec to ECI coordinates
            NRco0 = orbit.NR;
            NVco0 = orbit.NV;
            Nh_orbit = orbit.Nh;
            Nb2 = NRco0/norm(NRco0);
            Nb3 = -Nh_orbit/norm(Nh_orbit);
            Nb1 = cross(Nb2,Nb3);
            NQB = [Nb1 Nb2 Nb3];
            NRac0 = NQB*BRac0;
            NRbc0 = NQB*BRbc0;
            NRfc0 = NQB*BRfc0;
            NRgc0 = NQB*BRgc0;
            NRdc0 = NQB*BRdc0;
            NRec0 = NQB*BRec0;
            NVac0 = NQB*BVac0;
            NVbc0 = NQB*BVbc0;
            NVfc0 = NQB*BVfc0;
            NVgc0 = NQB*BVgc0;
            NVdc0 = NQB*BVdc0;
            NVec0 = NQB*BVec0;
            obj.NRao0 = NRac0+NRco0;
            obj.NRbo0 = NRbc0+NRco0;
            obj.NRfo0 = NRfc0+NRco0;
            obj.NRgo0 = NRgc0+NRco0;
            obj.NRdo0 = NRdc0+NRco0;
            obj.NReo0 = NRec0+NRco0;
            obj.NVao0 = NVac0+NVco0;
            obj.NVbo0 = NVbc0+NVco0;
            obj.NVfo0 = NVfc0+NVco0;
            obj.NVgo0 = NVgc0+NVco0;
            obj.NVdo0 = NVdc0+NVco0;
            obj.NVeo0 = NVec0+NVco0;
        end
    end
end

