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
        function obj = state_vec_N(BX, orbit)
            %STATE_VEC_N Construct an instance of this class
            %   Detailed explanation goes here
            NRco0 = orbit.NR;
            NVco0 = orbit.NV;
            Nh_orbit = orbit.Nh;
            Nb2 = -NRco0/norm(NRco0);
            Nb3 = Nh_orbit/norm(Nh_orbit);
            Nb1 = cross(Nb2,Nb3);
            NQB = [Nb1 Nb2 Nb3];
            NRac0 = NQB*BX.BRac0;
            NRbc0 = NQB*BX.BRbc0;
            NRfc0 = NQB*BX.BRfc0;
            NRgc0 = NQB*BX.BRgc0;
            NRdc0 = NQB*BX.BRdc0;
            NRec0 = NQB*BX.BRec0;
            NVac0 = NQB*BX.BVac0;
            NVbc0 = NQB*BX.BVbc0;
            NVfc0 = NQB*BX.BVfc0;
            NVgc0 = NQB*BX.BVgc0;
            NVdc0 = NQB*BX.BVdc0;
            NVec0 = NQB*BX.BVec0;
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

