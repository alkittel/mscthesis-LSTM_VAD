function q=v_rotax2qr(a,t)
%V_ROTQR2AX converts a v_rotation axis and angle to the corresponding real quaternion
% Inputs:
%
%     A(3,1)   Vector in the direction of the v_rotation axis.
%     T        Rotation angle in radians
%
% Output: 
%
%     Q(4,1)   normalized real-valued quaternion
%
% In the quaternion representation of a v_rotation, and q(1) = cos(t/2) 
% where t is the angle of v_rotation and q(2:4)/sin(t/2) is a unit vector
% lying along the axis of v_rotation.
% A positive v_rotation about [0 0 1] takes the X axis towards the Y axis.
% 
%      Copyright (C) Mike Brookes 2007-2018
%      Version: $Id: v_rotax2qr.m 10865 2018-09-21 17:22:45Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/v_voicebox/v_voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if all(a==0)
    a=[1;0;0];
end
m=sqrt(a(:)'*a(:));
q=[cos(0.5*t); sin(0.5*t)*a(:)/m];
if ~nargout
    v_rotqr2ro(q); % plot a rotated cube
end
    