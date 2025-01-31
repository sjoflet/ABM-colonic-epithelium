classdef BT3 % create object

    properties
        rep_rate ; 
        last_rep ; 
        biomass ; 
        position ;  
        alive ; 
    end 
    methods 
        
        function BT3 = BT3(varargin)


            switch nargin
                case 0 
                    BT3.rep_rate = [] ; 
                    BT3.last_rep = [] ;  
                    BT3.biomass = [] ; 
                    BT3.position = [] ; 
                    BT3.alive = [] ; 
                   
                case 1
                    if (isa(varargin{1}, 'BT3'))
                        BT3 = varargin{1} ; 
                    else
                        error('Input argument is not BT3')

                    end 
                case 5 
                    BT3.rep_rate = varargin{1} ; 
                    BT3.last_rep = varargin{2} ; 
                    BT3.biomass = varargin{3} ; 
                    BT3.position = varargin{4} ;
                    BT3.alive = varargin{5} ; 
                    
                otherwise 
                    error('Invalid no. of input arguments')
            end 
        end
    end
end
