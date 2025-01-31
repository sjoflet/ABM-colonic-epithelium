classdef BT1 % create object
    properties 

        rep_rate ; 
        last_rep ; 
        position ; 
        biomass  ; 
        anchor ; 
        alive ; 

    end 
    methods 
        
        function BT1 = BT1(varargin)
            switch nargin
                case 0 
                    BT1.rep_rate = [] ; 
                    BT1.last_rep= [] ; 
                    BT1.position = [] ; 
                    BT1.biomass = [] ; 
                    BT1.anchor = [] ;
                    BT1.alive = [] ; 
                case 1
                    if (isa(varargin{1}, 'BT1'))
                        BT1 = varargin{1} ; 
                    else
                        error('Input argument is not BT1')

                    end 
                case 6 
                    BT1.rep_rate = varargin{1};
                    BT1.last_rep= varargin{2} ; 
                    BT1.position = varargin{3}; 
                    BT1.biomass = varargin{4}; 
                    BT1.anchor = varargin{5};
                    BT1.alive = varargin{6}; 

                otherwise 
                    error('Invalid no. of input arguments')
            end 
        end
    end
end 