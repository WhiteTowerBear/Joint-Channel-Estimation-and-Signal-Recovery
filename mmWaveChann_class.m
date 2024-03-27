classdef mmWaveChann_class < handle 
    %    % mmWaveChann_class                                               
    % The class mmWaveChann_class will generate channel coefficients for
    % millimeter wave signals. Different models can be used. 
    % 
    % Saleh_Valenzuel Model (default)
    %   please read "Spatially Sparse Precoding in Millimeter Wave MIMO
    %   Systems" by Omar El Ayach for more information.
    %
    % A user can change any desired parameter in the constructor. 
    %   mmWaveChann_class(PropertyName0,PrpertyValue0,...
    %                     PropertyName1,PrpertyValue1,...
    %                     ...
    %                     PropertyNameN,PrpertyValueN)
    %
    %  Properties that can be set:
    %         ChannelType
    %         CarrierFrequency
    %         NumberOfMainPaths
    %         NumberOfSubPathsPerMainPath
    %         TxNumberOfAntennasInEachColumn
    %         TxNumberOfAntennasInEachRows
    %         TxAzimushMaxAngle
    %         TxAzimushMinAngle
    %         TxAzimushAngleSubPathStd
    %         TxElevationMaxAngle
    %         TxElevationMinAngle
    %         TxElevationAngleSubPathStd
    %         TxInterElementSpacing
    %         RxNumberOfAntennasInEachColumn
    %         RxNumberOfAntennasInEachRows
    %         RxAzimushMaxAngle
    %         RxAzimushMinAngle
    %         RxAzimushAngleSubPathStd
    %         RxElevationMaxAngle
    %         RxElevationMinAngle
    %         RxElevationAngleSubPathStd
    %         RxInterElementSpacing
    %         TxNumberOfAntennas
    %         RxNumberOfAntennas
    %         Lambda
    %         nTx
    %         nRx  
    % 
    %   see also: mmWaveChann_class.step, mmWaveChann_class.plot,
    %   mmWaveChann_class.TC001
    %  
    %  Examples 01:
    %     % construct an object with desired parameters
    %     c = mmWaveChann_class('TxNumberOfAntennasInEachRows'  ,4,...)
    %                           'TxNumberOfAntennasInEachColumn',4,...
    %                           'RxNumberOfAntennasInEachRows'  ,3,...
    %                           'RxNumberOfAntennasInEachColumn',3);
    %     % get channel coefficients
    %     H1 = c.step();    % generate channel coeff 1
    %     H2 = c.step();    % generate channel coeff 2 independent of 1
    %
    %     % run a test-case to see how random channel coeffs are constructed.
    %     c.TC001();
    %
    %     % run a test-case to checks the Frobenius norm of the channel.
    %     c.TC002();
    %
    %     % get help for a property
    %     help c.ChannelType
    %
    %
    %----------------------------------------------------------------------
    % Author    : Hamid Ramezani
    % Date      : 16-Oct-2017
    % email     : hamid.ramezani@gmail.com
    % Version   : 1.00 first time creation
    % MATLAB    : R2015b
    %----------------------------------------------------------------------
    
    properties (SetAccess = private) % Main properties --------------------
        
        ChannelType                     = 'Saleh-Valenzuel' % channel model: 'Saleh-Valenzuel'
        CarrierFrequency                = 26e9;    % carrier frequency in Hz
                
        NumberOfMainPaths               = 3;       % number of main paths (clusters)
        NumberOfSubPathsPerMainPath     = 8;       % number of sub-paths per main paths

        TxNumberOfAntennasInEachColumn  = 4        % number of antenna row     in a uniformly planar array Tx
        TxNumberOfAntennasInEachRows    = 4        % number of antenna columns in a uniformly planar array Tx
        TxAzimushMaxAngle               =  pi;     % for random angle generation
        TxAzimushMinAngle               = -pi;     % for random angle generation
        TxAzimushAngleSubPathStd        =  pi/128; % for each cluster
        TxElevationMaxAngle             =  pi/2;   % for random angle generation
        TxElevationMinAngle             = -pi/2;   % for random angle generation
        TxElevationAngleSubPathStd      =  pi/128; % for each cluster
        TxInterElementSpacing           =  3/520;  % space between two adjacent antennas (now lambda/2)
        
        RxNumberOfAntennasInEachColumn   = 4        % number of antenna row     in a uniformly planar array Rx     
        RxNumberOfAntennasInEachRows     = 1        % number of antenna columns in a uniformly planar array Rx 
        RxAzimushMaxAngle                =  pi;     % maximum value that an azimuth angle can get
        RxAzimushMinAngle                = -pi;     % minimum value that an azimuth angle can get
        RxAzimushAngleSubPathStd         =  pi/32;  % standard deviation of the scattered sub-path
        RxElevationMaxAngle              =  pi/2;   % maximum value that an elevation angle can get
        RxElevationMinAngle              = -pi/2;   % minimum value that an elevation angle can get
        RxElevationAngleSubPathStd       =  pi/32;  % standard deviation of the scattered sub-path
        RxInterElementSpacing           =   3/520;  % space between two adjacent antennas
                
    end
    properties (Dependent) % Dependent parameters -------------------------
        TxNumberOfAntennas  % number of transmit antennas      
        RxNumberOfAntennas  % number of receive  antennas
        Lambda              % wave length
        nTx                 % number of transmit antennas same as TxNumberOfAntennas
        nRx                 % number of receive  antennas same as RxNumberOfAntennas
    end
    methods % constructure ------------------------------------------------
        function c = mmWaveChann_class(varargin)%--------------------------
            % set variables
            for i = 1 : 2 : length(varargin)
                if isa(c.(varargin{i}),class(varargin{i+1}))
                    c.(varargin{i}) = varargin{i+1};
                else
                    error('Value for property "%s" has to be "%s"',...
                        varargin{i},class(c.(varargin{i})))
                end
            end            
        end
    end   
    methods % main functions ----------------------------------------------
        function [H,P] = step(c)%------------------------------------------
            % [H,P] = ().step();
            %  H : generated coefficients of the channel
            %  P : A struct which holds all the information about the
            %      subpaths, i.e., the AoA (azimuth and elevation), the  
            %      AoD (azimuth and elevation) and the gain of each path.
            %
            % see also mmWaveChann_class
            
            %#codegen
            % generate random angles for AoD (uniform distribution)
            txPhi   = rand(c.NumberOfMainPaths,1)*...
                (c.TxAzimushMaxAngle   - c.TxAzimushMinAngle)  +c.TxAzimushMinAngle;
            txTheta = rand(c.NumberOfMainPaths,1)*...
                (c.TxElevationMaxAngle - c.TxElevationMinAngle)+c.TxElevationMinAngle;

            rxPhi   = rand(c.NumberOfMainPaths,1)*...
                (c.RxAzimushMaxAngle   - c.RxAzimushMinAngle)  +c.RxAzimushMinAngle;
            rxTheta = rand(c.NumberOfMainPaths,1)*...
                (c.RxElevationMaxAngle - c.RxElevationMinAngle)+c.RxElevationMinAngle;

            gamma = sqrt((c.TxNumberOfAntennas*c.RxNumberOfAntennas) ./ ...
                         (c.NumberOfMainPaths *c.NumberOfSubPathsPerMainPath));                         

            % generate random complex numbers
            alpha = sqrt(1/2).* gamma * ...
                           (      randn(c.NumberOfMainPaths,c.NumberOfSubPathsPerMainPath) ...
                             + 1i*randn(c.NumberOfMainPaths,c.NumberOfSubPathsPerMainPath) );
            % coefficients of the subpath 
            txPhiSubPath   = zeros(c.NumberOfMainPaths,c.NumberOfSubPathsPerMainPath);
            txThetaSubPath = zeros(c.NumberOfMainPaths,c.NumberOfSubPathsPerMainPath);
            rxPhiSubPath   = zeros(c.NumberOfMainPaths,c.NumberOfSubPathsPerMainPath);
            rxThetaSubPath = zeros(c.NumberOfMainPaths,c.NumberOfSubPathsPerMainPath);
            for n = 1 : c.NumberOfMainPaths
                % generate random azimuf and elevation angles
                txPhiSubPath(n,:) = txPhi(n) + ...
                    randn(c.NumberOfSubPathsPerMainPath,1).*c.TxAzimushAngleSubPathStd;
                txThetaSubPath(n,:) = txTheta(n) + ...
                    randn(c.NumberOfSubPathsPerMainPath,1).*c.TxElevationAngleSubPathStd;
                rxPhiSubPath(n,:) = rxPhi(n) + ...
                    randn(c.NumberOfSubPathsPerMainPath,1).*c.RxAzimushAngleSubPathStd;
                rxThetaSubPath(n,:) = rxTheta(n) + ...
                    randn(c.NumberOfSubPathsPerMainPath,1).*c.RxElevationAngleSubPathStd;
            end            
            % Initialize H with zeros
            H = complex(zeros(c.RxNumberOfAntennas,c.TxNumberOfAntennas));
            
            if strcmpi(c.ChannelType,'Saleh-Valenzuel')
                for n = 1 : c.NumberOfMainPaths
                    % generate random azimuth and elevation angles
                    for m = 1 : c.NumberOfSubPathsPerMainPath
                        aTx = c.txBuildArrayReponses(txPhiSubPath(n,m),txThetaSubPath(n,m));
                        aRx = c.rxBuildArrayReponses(rxPhiSubPath(n,m),rxThetaSubPath(n,m));
                        H = H + ...
                           alpha(n,m) * aRx * aTx';
                    end
                end
            else
                error('Such channel type is not supported.')
            end
            if nargout == 2
                P.txPhiSubPath   = txPhiSubPath;
                P.txThetaSubPath = txThetaSubPath;
                P.rxPhiSubPath   = rxPhiSubPath;
                P.rxThetaSubPath = rxThetaSubPath;
                P.alpha          = alpha;
                P.H              = H;
            end
        end
        function aTx = txBuildArrayReponses(c,phi,theta)%------------------
            % for planar case which is general
            [n,m] = meshgrid(0:(c.TxNumberOfAntennasInEachRows-1),...
                             0:(c.TxNumberOfAntennasInEachColumn-1));                         
            aTx   = 1./sqrt(c.TxNumberOfAntennas) .* ...
                exp(2*pi*c.TxInterElementSpacing/c.Lambda*1i*(m(:).*sin(phi)*sin(theta)+n(:).*cos(theta)));
        end
        function aRx = rxBuildArrayReponses(c,phi,theta)%------------------
            % for planar case which is general
            [n,m] = meshgrid(0:(c.RxNumberOfAntennasInEachRows-1),...
                             0:(c.RxNumberOfAntennasInEachColumn-1));
            aRx   = 1./sqrt(c.RxNumberOfAntennas) .* ...
                exp(2*pi*c.RxInterElementSpacing/c.Lambda*1i*(m(:).*sin(phi)*sin(theta)+n(:).*cos(theta)));
        end
    end
    methods % for dependent parameters ------------------------------------
        function value = get.TxNumberOfAntennas(c)%------------------------
            value = c.TxNumberOfAntennasInEachRows * ...
                    c.TxNumberOfAntennasInEachColumn;
        end
        function value = get.RxNumberOfAntennas(c)%------------------------
            value = c.RxNumberOfAntennasInEachRows * ...
                    c.RxNumberOfAntennasInEachColumn;
        end   
        function value = get.Lambda(c)%------------------------------------
            value = 3e8/c.CarrierFrequency;
        end
    end
    methods % related to plots --------------------------------------------
        function plot(c,P)%------------------------------------------------
            % ().plot(P)
            %  Plots channel properties.
            %
            % P is a struct holding all channel information, generated by 
            % [H,P] = ().step();
            %
            %
            figure(1)
            clf
            subplot(2,2,1)
            for n = 1 : c.NumberOfMainPaths
                % for transmitter
                plot(cos(P.txPhiSubPath(n,:)),sin(P.txPhiSubPath(n,:)),'o')
                hold on 
                plot(cos(P.rxPhiSubPath(n,:)),sin(P.rxPhiSubPath(n,:)),'x')                
            end
            x = linspace(-1,1,100);
            y = sqrt(1-x.^2);
            plot([x,x(end:-1:1)],[y,-y(end:-1:1)],'k:')
            axis square
            xlabel('cos(\phi)')
            ylabel('sin(\phi)')
            grid on 
            legend('tx \phi','rx \phi','location','best')
            title('azimuth angles of each path from tx to rx')
            
            subplot(2,2,2)
            for n = 1 : c.NumberOfMainPaths
                % for transmitter
                plot(cos(P.txThetaSubPath(n,:)),sin(P.txThetaSubPath(n,:)),'o')
                hold on 
                plot(cos(P.rxThetaSubPath(n,:)),sin(P.rxThetaSubPath(n,:)),'x')                
            end
            x = linspace(-1,1,100);
            y = sqrt(1-x.^2);
            plot([x,x(end:-1:1)],[y,-y(end:-1:1)],'k:')
            axis square
            xlabel('cos(\theta)')
            ylabel('sin(\theta)')
            grid on 
            legend('tx \theta','rx \theta','location','best')
            title('elevation angles of each path from tx to rx')

            subplot(2,2,3)
            bar(abs(P.alpha))
            axis square
            xlabel('path index')
            ylabel('path amplitude')
            grid on 
            title('amplitude of each path-subpath')
            
            subplot(2,2,4)
            S = svd(P.H);
            plot(S)
            xlabel('singular index')
            ylabel('singular value')    
            grid on
            axis square
            title('singula values of the channel')
            
        end
    end
    methods % related to test-cases ---------------------------------------
        function TC001(c)% ------------------------------------------------
            % A test case to show how the channel is built from AoAs, 
            % AoDs and gains
            
            % set random seed
            rng(0)
            % generate a channel 
            [~,P] = c.step;
            % plot channel 
            c.plot(P);
        end
        function TC002(c)% ------------------------------------------------
            % checks the Frobenius norm of the channel
            rng(0)
            N = 1e4;
            O = zeros(N,1);
            for n = 1 : N
                H    = c.step();
                O(n) = norm(H(:))^2;
            end
            % average the norm
            Average_norm = mean(O);
            
            % display results
            fprintf('Simulation E{|H(:)|^2_F} = %3.2f \n',Average_norm)
            fprintf('Required   E{|H(:)|^2_F} = %3.2f \n',c.TxNumberOfAntennas*c.RxNumberOfAntennas)
        end
    end
end