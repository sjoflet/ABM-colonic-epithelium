classdef Goblet % create object
properties 

        muc2_prodrate ; 
        muc2_lastProd ; 
        position ; 
        butyrateTot ; 
        alive ; 
        

    end 
    methods 
        
        function gob = Goblet(varargin)


            switch nargin
                case 0 
                    gob.muc2_prodrate = [] ; 
                    gob.position = [] ;  
                    gob.alive = [] ; 
                case 1
                    if (isa(varargin{1}, 'Bacteria1'))
                        gob = varargin{1} ; 
                    else
                        error('Input argument is not Bacteria 1')

                    end 
                case 5
                    gob.muc2_prodrate = varargin{1} ;
                    gob.muc2_lastProd = varargin{2} ; 
                    gob.position = varargin{3} ;
                    gob.butyrateTot = varargin{4} ; 
                    gob.alive = varargin{5} ; 

                otherwise 
                    error('Invalid no. of input arguments')
            end 
        end
    end



end 