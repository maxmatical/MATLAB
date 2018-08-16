function [lsvec,rsvec,sig] = get_dom_triple(A,tolr,tolc, cdict, files)
%% [lsvec,rsvec,sig] = get_dom_triple(A,tolr,tolc, cdict, files)
%% Find dominant singular triple of a submatrix of A.  Input arguments:
%%  tolr= angle tolerance for dropping a row of A from the submatrix
%%  tolc = angle tolerance for allowing a column of A to enter the
%%  submatrix
%% cdict,files - unused arguments; available for debugging purposes.

nr = norm(A,'fro');
[m,n] = size(A);
[scrap,ipos] = max(sum(A));
rsvec_sp = [1];
sig = 0;

%% Initialize the columns we are keeping to just one (ipos) in cmask.
%% Initialize the rows we are keeping to all of them in rmask.

cmask = logical(zeros(n,1));
cmask(ipos) = 1;
rmask = logical(ones(m,1));
Atran = A';
stabcount = 0;


cutoff_decr = 0.05;
cutoff_incr = 0.02;

current_rcutoff = 0.01;
current_ccutoff = 0.95;

while 1    
  row_norms = sqrt(sum(Atran(cmask,:).^2,1))';
  old_rmask = rmask;    
  
  % Update the left singular vector.
  
  lsvec0 = A(:,cmask) * rsvec_sp;    
  
  % Compute the row angles.
  angs = abs(lsvec0) ./ (max(row_norms,eps) * norm(rsvec_sp));
  toss_flag = angs < current_rcutoff;
  toss_flag(row_norms == 0) = 0;
  rmask(toss_flag) = 0;  
  
  
  lsvec_sp = lsvec0(rmask);
  lsvec_sp = lsvec_sp / norm(lsvec_sp); 
  
  col_norms = sqrt(sum(A(rmask,:).^2,1))';
  old_cmask = cmask;
  
  %% Update the right singular vector
  rsvec0 = Atran(:,rmask) * lsvec_sp;  
  
  % Compute the column angles
  angs2 = abs(rsvec0) ./ (max(col_norms,eps) * norm(lsvec_sp));
  
  keep_flag2 = angs2 >= current_ccutoff & col_norms > 0;
  cmask(keep_flag2) = 1;
  rsvec_sp = rsvec0(cmask);
  
   
  
  oldsig = sig;
  sig = norm(rsvec_sp) / norm(lsvec_sp);
  
  rsvec_sp = rsvec_sp / norm(rsvec_sp);
  
  
  %% Termination criterion: cmask unchanged, rmask unchanged,
  %% the singular value is converging by a relative tolerance of 1e-6,
  %% three successive iterations with these conditions holding.
  
  if abs(oldsig - sig) < 1e-6 * nr && all(old_rmask == rmask) && ...
      all(old_cmask == cmask) && current_rcutoff >= tolr && current_ccutoff <= tolc
    if stabcount == 3
      break
    else
      stabcount = stabcount + 1;
    end
  else
    stabcount = 0;
  end
  if current_rcutoff < tolr
    current_rcutoff = current_rcutoff + cutoff_incr;
  end
  if current_ccutoff > tolc
    current_ccutoff = current_ccutoff - cutoff_decr;
  end  
end

lsvec = zeros(m,1);
lsvec(rmask) = lsvec_sp;

rsvec = zeros(n,1);
rsvec(cmask) = rsvec_sp;
disp(sprintf(...
  'dom_triple initialized with col %d; final # of cols = %d # of rows = %d', ...
  ipos, nnz(rsvec_sp), nnz(lsvec_sp)))