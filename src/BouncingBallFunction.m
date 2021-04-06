function BouncingBallFunction
 
%--------------------------------------------------------------------------
%BouncingBallFunction
%Creating Simple Animation in MATLAB
%MATLAB Undercover
%zerocrossraptor.wordpress.com
%--------------------------------------------------------------------------
%This function m-file creates an endless animation of bouncing ball.
%--------------------------------------------------------------------------
 
%CodeStart-----------------------------------------------------------------
%Declaring ball's initial condition
    initpos=50;     %Ball's initial vertical position
    initvel=0;      %Ball's initial vertical velocity
%Declaring environmental variable
    r_ball=5;       %Ball's radius
    gravity=10;     %Gravity's acceleration
    c_bounce=1;     %Bouncing's coefficient of elasticity
%Declaring animation timestep
    dt=0.0125;      %Animation timestep
%Initiating figure, axes, and objects for animation
    fig=figure('DeleteFcn',@closefigurefcn);
    axs=axes('Parent',fig);
    ball=rectangle('Position',[-r_ball,initpos,r_ball,r_ball],...
                   'Curvature',[1,1],...
                   'FaceColor','b',...
                   'Parent',axs);
    line([-5*r_ball,5*r_ball],...
         [0,0],...
         'Parent',axs);
%Executing animation
    pos=initpos-r_ball;             %Ball's current vertical position
    vel=initvel;                    %Ball's current vertical velocity
    play=true;                      %Current animation status
    while play
        %Declaring time counter
        t_loopstart = tic();
        %Updating ball's condition
        pos=pos+(vel*dt);           %Ball's current vertical position
        vel=vel-(gravity*dt);       %Ball's current vertical velocity
        if pos<0
            vel=-vel*c_bounce;      %Ball's current vertical velocity
        end
        %Updating ball
        set(ball,'Position',[-r_ball,pos,r_ball,r_ball]);
        %Preserving axes
        axis(axs,[-5*r_ball,5*r_ball,0,initpos+2*r_ball]);
        axis(axs,'equal');
        axis(axs,'off');
        %Pausing animation
        el_time=toc(t_loopstart);
        disp(['Elapse time : ',num2str(el_time),' seconds']);
        pause(dt-el_time);
    end
%Declaring callback function
    function closefigurefcn(source,event)
        play=false;
        pause(1);
    end
%CodeEnd-------------------------------------------------------------------
 
end