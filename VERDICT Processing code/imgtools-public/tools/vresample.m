function B = vresample(vV,mVi,mQi) 
% VRESAMPLE Resample volume using DICOM-style geom structures
%
% B = vresample(vV,mV,mQ) 
% B = vresample(vV,Vgeom,Qgeom) 
%
% Refactor of previous dreslice function
%
% Example usage:
%   dinfoV = dmfparse(dselect) ;
%   dinfoQ = dmfparse(dselect) ;
%
%   [vV, mV] = d2mat(dinfoV,{'slice'},'op','fp');
%   [vQ, mQ] = d2mat(dinfoQ,{'slice'},'op','fp');
%
%   Optional: to change the pixel spacing in the output:
%   mQ.geom = geom_change_ps(mQ.geom, mV.geom(1).PixelSpacing_HW) 
%
%   B = vresample(vV,mV,mQ) ;
%
% THIS SOFTWARE IS PROVIDED  "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% Copyright 2020, David Atkinson
% D.Atkinson@ucl.ac.uk
%
% See also SET_INTERPOLANT DRESLICE GEOM_CHANGE_PS  dicom2intrinsic

if ~isfield(mQi,'geom'), mQ.geom = mQi ; else, mQ = mQi ;  end
if ~isfield(mVi,'geom'), mV.geom = mVi ; else, mV = mVi ; end

nQslice = length(mQ.geom) ;

AV2i = dicom2intrinsic(mV.geom) ; % maps V's LPH to reference intrinsic
AQ2i = dicom2intrinsic(mQ.geom) ;

tformAQD2i = affine3d(AQ2i) ; % tform class (not struct) 
tformAVD2i = affine3d(AV2i) ;

% get Q's own intrinsic coordimates (these are not aligned with Vs)
% Using meshgrid (not ndgrid).
% QI1 will come out with the expected shape - which flows through to tmap_B
% and B.
% QI1 will be the 'x' coordinate values
[QI1, QI2, QI3] = meshgrid(1:mQ.geom(1).Width, 1:mQ.geom(1).Height, 1:nQslice) ;

% Get the DICOM LPH coordinates of Q
[QD1, QD2, QD3] = transformPointsInverse(tformAQD2i, QI1, QI2, QI3) ;

% Convert the DICOM coordinates of Q to the intrinsic of V (ready for
% interpolation within V using tformarray)
[QVI1, QVI2, QVI3] = transformPointsForward(tformAVD2i, QD1, QD2, QD3) ;

% tmap_B are the positions of the points in Q expressed in terms of the
% intrinsic coordinates of V
tmap_B = cat(4,QVI1, QVI2, QVI3) ;

interpolant = set_interpolant(mV.geom, mQ.geom) ;

R = makeresampler(interpolant,'fill') ;

% The [2 1 3] 'swap' I think is because the order in tmap_B is X,Y,Z so it
% needs to know that the mapping to vV requires a swap??
B = tformarray(vV, [], R, [2 1 3],[1 2 3],[],tmap_B,0);
