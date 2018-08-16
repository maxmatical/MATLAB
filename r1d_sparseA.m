function [P,Q] = r1d_sparseA(A,termlist,files,k,tolr,tolc)
% Run [P,Q] = r1d_sparseA(A,[],[],k,tolr,tolc)
%k: the number of bases 
%%  tolr= angle tolerance for dropping a row of A from the submatrix
%%  tolc = angle tolerance for allowing a column of A to enter the
[m,n] = size(A);
P = sparse([],[],[],m,k);
Q = sparse([],[],[],n,k);
for j = 1 : k 
  
  [u,v,sig] = get_dom_triple(A,tolr,tolc,termlist,files);
  %% The remainder of this loop executes the statements
  %% A = max(A-u*sig*v', zeros(m,n))
  %% P(:,j)=u*sig;
  %% Q(:,j)=v;
  %% with lots of sparse matrix tricks for efficiency.
  
  [A_is,A_js,A_es] = find(A);
  nnzA = length(A_is);
  [usparse_i, scrap, usparse_e] = find(u);
  nnzu = length(usparse_e);
  [vsparse_i, scrap, vsparse_e] = find(v);
  nnzv = length(vsparse_e);
  P(:,j) = sparse(usparse_i, ones(nnzu,1), usparse_e*sig,m,1);
  Q(:,j) = sparse(vsparse_i, ones(nnzv,1), vsparse_e,n,1);
  
  %% Find the common entries between A and u*v'.  Two ways to do this
  %% depending on whether A is sparser or u*v' is sparser.
  %% In either case, output the following vectors, which all have the
  %% same length (one vector entry per common matrix entry):
  %%  common_is: the row index of the common entry
  %%  common_js: the col index of the common entry
  %%  common_uve: the entry of u*v' corresponding to the common entry.
  %%  common_Ae: the entry of A corresponding to the common entry.
  
  if nnzu * nnzv < nnzA
    uv_is = kron(ones(nnzv,1),usparse_i);
    uv_js = kron(vsparse_i, ones(nnzu,1));
    uv_es = kron(vsparse_e,usparse_e);
    unified_Aind = (A_js - 1) * m + A_is;
    unified_uvind = (uv_js - 1) * m + uv_is;
    pos1 = ismembc2(unified_uvind, unified_Aind);
    keep_entries = find(pos1 > 0);
    common_is = uv_is(keep_entries);
    common_js = uv_js(keep_entries);
    common_uve = uv_es(keep_entries);
    common_Ae = A_es(pos1(keep_entries));
  else
    rowmatch = ismembc2(A_is, usparse_i);
    colmatch = ismembc2(A_js, vsparse_i);
    bothmatch = find((rowmatch > 0) & (colmatch > 0));
    common_is = A_is(bothmatch);
    common_js = A_js(bothmatch);
    common_uve = u(A_is(bothmatch)) .* v(A_js(bothmatch));
    common_Ae = A_es(bothmatch);
  end  
  update_e = min(sig*common_uve, common_Ae);
  A = A - sparse(common_is, common_js, update_e, m, n);
end