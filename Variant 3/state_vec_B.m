classdef state_vec_B
    %STATE_VEC_BCOORDS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        BWc0
        BRac0
        BRbc0
        BRfc0
        BRgc0
        BRdc0
        BRec0
        BVac0
        BVbc0
        BVfc0
        BVgc0
        BVdc0
        BVec0
    end
    
    methods
        function obj = state_vec_B(p)
            %STATE_VEC_BCOORDS Construct an instance of this class
            %   Detailed explanation goes here
            b1 = [1;0;0]; % basis vectors
            b2 = [0;1;0];
            b3 = [0;0;1];
            obj.BWc0 = p.w*b3;
            obj.BRac0 = p.re*b1;
            obj.BRbc0 = -p.re*b1;
            obj.BRfc0 = p.re*b2;
            obj.BRgc0 = -p.re*b2;
            obj.BRdc0 = p.de*b3;
            obj.BRec0 = -p.de*b3;
            obj.BVac0 = cross(obj.BWc0, obj.BRac0);
            obj.BVbc0 = cross(obj.BWc0, obj.BRbc0);
            obj.BVfc0 = cross(obj.BWc0, obj.BRfc0);
            obj.BVgc0 = cross(obj.BWc0, obj.BRgc0);
            obj.BVdc0 = zeros(3, 1);
            obj.BVec0 = zeros(3, 1);
        end
    end
end

