function [t,VAR,Output] = C172stability
%===========================================================================
% File: C172stability.m created Dec 20 2022 by MotionGenesis 6.1.
% Portions copyright (c) 2009-2021 Motion Genesis LLC.  Rights reserved.
% MotionGenesis Student Licensee: Maxime Couture. (until September 2026).
% Paid-up MotionGenesis Student licensees are granted the right
% to distribute this code for legal student-academic (non-professional) purposes only,
% provided this copyright notice appears in all copies and distributions.
%===========================================================================
% The software is provided "as is", without warranty of any kind, express or    
% implied, including but not limited to the warranties of merchantability or    
% fitness for a particular purpose. In no event shall the authors, contributors,
% or copyright holders be liable for any claim, damages or other liability,     
% whether in an action of contract, tort, or otherwise, arising from, out of, or
% in connection with the software or the use or other dealings in the software. 
%===========================================================================
eventDetectedByIntegratorTerminate1OrContinue0 = [];
IPlanezz=0; qaDDt=0; xDDt=0; yDDt=0; alpha=0; alpha_stab=0; Cd_stab=0; Cd_wing=0; Cl_stab=0; Cl_wing=0; Drag_stab=0; Drag_wing=0;
Lift_stab=0; Lift_wing=0;

Plane= 'C172'

%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
%deflection                      = -5;                      % degrees             Constant
Density                         =  1.225;                  % kg/m^3              Constant
efficiency                      =  0.735;                  % noUnits             Constant
g                               =  9.80665;                % m/s^2               Constant

power_a                         =  116.6;                  % kiloWatt            Constant
r                               =  2;                      % m                   Constant

if Plane == 'C172':

   S_stab                          =  3.96;                   % m^2                 Constant
   S_wing                          =  16.2;                   % m^2                 Constant
   x_stab                          = -4.24;                   % m                   Constant
   x_wing                          =  0.25;                   % m                   Constant
   y_stab                          =  0.5;                    % m                   Constant
   y_wing                          =  1;                      % m                   Constant
   mPlane                          =  779;                    % kg                  Constant

qa                              =  0;                      % degrees             Initial Value
x                               =  0;                      % m                   Initial Value
y                               =  1000;                   % m                   Initial Value
qaDt                            =  0;                      % UNITS               Initial Value
xDt                             =  100;                    % m/s                 Initial Value
yDt                             =  0;                      % m/s                 Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  10;                     % second              Final Time
tStep                           =  0.02;                   % sec                 Integration Step
printIntScreen                  =  1;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  1;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------


for deflection = 10:-10

   % Unit conversions
   DEGtoRAD = pi / 180.0;
   RADtoDEG = 180.0 / pi;
   deflection = deflection * DEGtoRAD;
   power_a = power_a * 1000;
   qa = qa * DEGtoRAD;

   % Evaluate constants
   IPlanezz = 0.5*mPlane*r^2;


   VAR = SetMatrixFromNamedQuantities;
   [t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
   OutputToScreenOrFile( [], 0, 0 );   % Close output files.
   if( printIntFile ~= 0 ),  PlotOutputFiles;  end


   %===========================================================================
   function sys = mdlDerivatives( t, VAR, uSimulink )
   %===========================================================================
   SetNamedQuantitiesFromMatrix( VAR );
   alpha = sign0IsPositive1(sin(qa)*xDt-cos(qa)*yDt)*acos((sin(qa)*yDt+cos(qa)*xDt)/sqrt(xDt^2+yDt^2));
   alpha_stab = deflection + alpha;
   Cd_stab = 2*sin(alpha_stab)^2;
   Drag_stab = 0.5*Density*S_stab*(xDt^2+yDt^2)*Cd_stab;
   Cd_wing = 2*sin(alpha)^2;
   Drag_wing = 0.5*Density*S_wing*(xDt^2+yDt^2)*Cd_wing;
   Cl_stab = 2*sin(alpha_stab)*cos(alpha_stab)*IsPositive(50-alpha);
   Lift_stab = 0.5*Density*S_stab*(xDt^2+yDt^2)*Cl_stab;
   Cl_wing = 2*sin(alpha)*cos(alpha)*IsPositive(50-alpha);
   Lift_wing = 0.5*Density*S_wing*(xDt^2+yDt^2)*Cl_wing;
   xDDt = 0.015625*(efficiency*power_a*cos(qa)-64*(xDt*Drag_stab+xDt*Drag_wing+yDt*Lift_stab+yDt*Lift_wing)/sqrt(xDt^2+yDt^2))/mPlane;
   yDDt = 0.015625*(efficiency*power_a*sin(qa)+64*(xDt*Lift_stab+xDt*Lift_wing-yDt*Drag_stab-yDt*Drag_wing)/sqrt(xDt^2+yDt^2))/mPlane - g;
   qaDDt = (x_stab*(sin(qa)*(xDt*Drag_stab+yDt*Lift_stab)+cos(qa)*(xDt*Lift_stab-yDt*Drag_stab))+x_wing*(sin(qa)*(xDt*Drag_wing+yDt*  ...
   Lift_wing)+cos(qa)*(xDt*Lift_wing-yDt*Drag_wing))+y_stab*(cos(qa)*(xDt*Drag_stab+yDt*Lift_stab)-sin(qa)*(xDt*Lift_stab-yDt*  ...
   Drag_stab))+y_wing*(cos(qa)*(xDt*Drag_wing+yDt*Lift_wing)-sin(qa)*(xDt*Lift_wing-yDt*Drag_wing)))/(IPlanezz*sqrt(xDt^2+yDt^2));

   sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
   end



   %===========================================================================
   function VAR = SetMatrixFromNamedQuantities
   %===========================================================================
   VAR = zeros( 1, 6 );
   VAR(1) = qa;
   VAR(2) = x;
   VAR(3) = y;
   VAR(4) = qaDt;
   VAR(5) = xDt;
   VAR(6) = yDt;
   end


   %===========================================================================
   function SetNamedQuantitiesFromMatrix( VAR )
   %===========================================================================
   qa = VAR(1);
   x = VAR(2);
   y = VAR(3);
   qaDt = VAR(4);
   xDt = VAR(5);
   yDt = VAR(6);
   end


   %===========================================================================
   function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
   %===========================================================================
   VARp = zeros( 1, 6 );
   VARp(1) = qaDt;
   VARp(2) = xDt;
   VARp(3) = yDt;
   VARp(4) = qaDDt;
   VARp(5) = xDDt;
   VARp(6) = yDDt;
   end



   %===========================================================================
   function Output = mdlOutputs( t, VAR, uSimulink )
   %===========================================================================
   Output = zeros( 1, 11 );
   Output(1) = t;
   Output(2) = x;
   Output(3) = y;
   Output(4) = qa*RADtoDEG;
   Output(5) = alpha*RADtoDEG;
   Output(6) = 0.0*RADtoDEG;

   Output(7) = x;
   Output(8) = y;

   Output(9) = t;
   Output(10) = alpha*RADtoDEG;
   Output(11) = alpha_stab*RADtoDEG;
   end


   %===========================================================================
   function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
   %===========================================================================
   persistent FileIdentifier hasHeaderInformationBeenWritten;

   if( isempty(Output) ),
      if( ~isempty(FileIdentifier) ),
         for( i = 1 : 3 ),  fclose( FileIdentifier(i) );  end
         clear FileIdentifier;
         fprintf( 1, '\n Output is in the files C172stability.i  (i=1,2,3)\n\n' );
      end
      clear hasHeaderInformationBeenWritten;
      return;
   end

   if( isempty(hasHeaderInformationBeenWritten) ),
      if( shouldPrintToScreen ),
         fprintf( 1,                '%%       t              x              y             qa            alpha         stab_w\n' );
         fprintf( 1,                '%%     (sec)           (m)            (m)         (degrees)      (degrees)    (degrees/sec)\n\n' );
      end
      if( shouldPrintToFile && isempty(FileIdentifier) ),
         FileIdentifier = zeros( 1, 3 );
         FileIdentifier(1) = fopen('C172stability.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file C172stability.1'); end
         fprintf(FileIdentifier(1), '%% FILE: C172stability.1\n%%\n' );
         fprintf(FileIdentifier(1), '%%       t              x              y             qa            alpha         stab_w\n' );
         fprintf(FileIdentifier(1), '%%     (sec)           (m)            (m)         (degrees)      (degrees)    (degrees/sec)\n\n' );
         FileIdentifier(2) = fopen('C172stability.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file C172stability.2'); end
         fprintf(FileIdentifier(2), '%% FILE: C172stability.2\n%%\n' );
         fprintf(FileIdentifier(2), '%%       x              y\n' );
         fprintf(FileIdentifier(2), '%%      (m)            (m)\n\n' );
         FileIdentifier(3) = fopen('C172stability.3', 'wt');   if( FileIdentifier(3) == -1 ), error('Error: unable to open file C172stability.3'); end
         fprintf(FileIdentifier(3), '%% FILE: C172stability.3\n%%\n' );
         fprintf(FileIdentifier(3), '%%       t            alpha       alpha_stab\n' );
         fprintf(FileIdentifier(3), '%%     (sec)        (degrees)      (degrees)\n\n' );
      end
      hasHeaderInformationBeenWritten = 1;
   end

   if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:6) );  end
   if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:6) );  end
   if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(7:8) );  end
   if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(3), Output(9:11) );  end
   end


   %===========================================================================
   function WriteNumericalData( fileIdentifier, Output )
   %===========================================================================
   numberOfOutputQuantities = length( Output );
   if( numberOfOutputQuantities > 0 ),
      for( i = 1 : numberOfOutputQuantities ),
         fprintf( fileIdentifier, ' %- 14.6E', Output(i) );
      end
      fprintf( fileIdentifier, '\n' );
   end
   end



   %===========================================================================
   function PlotOutputFiles
   %===========================================================================

   figure;
   data = load( 'C172stability.2' ); 
   plot( data(:,1),data(:,2),'-b', 'LineWidth',3 );
   legend( 'y (m)' );
   xlabel('x (m)');   ylabel('y (m)');   % title('Some plot title');
   clear data;

   figure;
   data = load( 'C172stability.3' ); 
   plot( data(:,1),data(:,2),'-b', data(:,1),data(:,3),'-.g', 'LineWidth',3 );
   legend( 'alpha (degrees)', 'alpha_stab (degrees)' );
   xlabel('t (sec)');   % ylabel('Some y-axis label');   title('Some plot title');
   clear data;
   end



   %===========================================================================
   function [functionsToEvaluateForEvent, eventTerminatesIntegration1Otherwise0ToContinue, eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1] = EventDetection( t, VAR, uSimulink )
   %===========================================================================
   % Detects when designated functions are zero or cross zero with positive or negative slope.
   % Step 1: Uncomment call to mdlDerivatives and mdlOutputs.
   % Step 2: Change functionsToEvaluateForEvent,                      e.g., change  []  to  [t - 5.67]  to stop at t = 5.67.
   % Step 3: Change eventTerminatesIntegration1Otherwise0ToContinue,  e.g., change  []  to  [1]  to stop integrating.
   % Step 4: Change eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1,  e.g., change  []  to  [1].
   % Step 5: Possibly modify function EventDetectedByIntegrator (if eventTerminatesIntegration1Otherwise0ToContinue is 0).
   %---------------------------------------------------------------------------
   % mdlDerivatives( t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
   % mdlOutputs(     t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
   functionsToEvaluateForEvent = [];
   eventTerminatesIntegration1Otherwise0ToContinue = [];
   eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1 = [];
   eventDetectedByIntegratorTerminate1OrContinue0 = eventTerminatesIntegration1Otherwise0ToContinue;
   end


   %===========================================================================
   function [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents )
   %===========================================================================
   isIntegrationFinished = eventDetectedByIntegratorTerminate1OrContinue0( nIndexOfEvents );
   if( ~isIntegrationFinished ),
      SetNamedQuantitiesFromMatrix( VAR );
   %  Put code here to modify how integration continues.
      VAR = SetMatrixFromNamedQuantities;
   end
   end



   %===========================================================================
   function [t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile )
   %===========================================================================
   OdeMatlabOptions = odeset( 'RelTol',relError, 'AbsTol',absError, 'MaxStep',tStep, 'Events',@EventDetection );
   t = tInitial;                 epsilonT = 0.001*tStep;                   tFinalMinusEpsilonT = tFinal - epsilonT;
   printCounterScreen = 0;       integrateForward = tFinal >= tInitial;    tAtEndOfIntegrationStep = t + tStep;
   printCounterFile   = 0;       isIntegrationFinished = 0;
   mdlDerivatives( t, VAR, 0 );
   while 1,
      if( (integrateForward && t >= tFinalMinusEpsilonT) || (~integrateForward && t <= tFinalMinusEpsilonT) ), isIntegrationFinished = 1;  end
      shouldPrintToScreen = printIntScreen && ( isIntegrationFinished || printCounterScreen <= 0.01 );
      shouldPrintToFile   = printIntFile   && ( isIntegrationFinished || printCounterFile   <= 0.01 );
      if( isIntegrationFinished || shouldPrintToScreen || shouldPrintToFile ),
         Output = mdlOutputs( t, VAR, 0 );
         OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile );
         if( isIntegrationFinished ), break;  end
         if( shouldPrintToScreen ), printCounterScreen = printIntScreen;  end
         if( shouldPrintToFile ),   printCounterFile   = printIntFile;    end
      end
      [TimeOdeArray, VarOdeArray, timeEventOccurredInIntegrationStep, nStatesArraysAtEvent, nIndexOfEvents] = ode45( @mdlDerivatives, [t tAtEndOfIntegrationStep], VAR, OdeMatlabOptions, 0 );
      if( isempty(timeEventOccurredInIntegrationStep) ),
         lastIndex = length( TimeOdeArray );
         t = TimeOdeArray( lastIndex );
         VAR = VarOdeArray( lastIndex, : );
         printCounterScreen = printCounterScreen - 1;
         printCounterFile   = printCounterFile   - 1;
         if( abs(tAtEndOfIntegrationStep - t) >= abs(epsilonT) ), warning('numerical integration failed'); break;  end
         tAtEndOfIntegrationStep = t + tStep;
         if( (integrateForward && tAtEndOfIntegrationStep > tFinal) || (~integrateForward && tAtEndOfIntegrationStep < tFinal) ) tAtEndOfIntegrationStep = tFinal;  end
      else
         t = timeEventOccurredInIntegrationStep( 1 );    % time  at firstEvent = 1 during this integration step.
         VAR = nStatesArraysAtEvent( 1, : );             % state at firstEvent = 1 during this integration step.
         printCounterScreen = 0;
         printCounterFile   = 0;
         [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents(1) );
      end
   end
   end



   %===========================================================================
   function retValue = sign0IsPositive1( x )
      if( x >= 0 ), retValue = 1;  else retValue = -1; end
   end
   function retValue = IsPositive( x )
      if( x > 0 ), retValue = 1;  else retValue = 0; end
   end

%======================================
end    % End of function C172stability
%======================================
