classdef BT2
properties 

        rep_rate ; 
        last_rep ; 
        position ; 
        biomass ; 
        alive ; 

    end 
    methods 
        
        function BT2 = BT2(varargin)


            switch nargin
                case 0 
                    BT2.rep_rate = [] ; 
                    BT2.last_rep= [] ; 
                    BT2.position = [] ; 
                    BT2.biomass = [] ;
                    BT2.alive = [ ]; 
                case 1
                    if (isa(varargin{1}, 'BT2'))
                        BT2 = varargin{1} ; 
                    else
                        error('Input argument is not E coli')

                    end 
                case 5 
                    BT2.rep_rate = varargin{1} ; 
                    BT2.last_rep= varargin{2} ; 
                    BT2.position = varargin{3} ; 
                    BT2.biomass = varargin{4} ; 
                    BT2.alive = varargin{5} ; 
                otherwise 
                    error('Invalid no. of input arguments')
            end 
        end
    end

end 
