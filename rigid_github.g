Sigma_gens := function(n)
	# Return a set of generators for Sylow 2-subgroup of Sym(2^)
    local dim, x, sgens, i, j;
    dim:=n; # the dimension n 
    sgens:=[]; # will contain generators for Sigma_n
    for i in [1..dim] do
	x:=();
        for j in [1..2^(i-1)] do
            x:=x*(j,j+2^(i-1));
        od;
        Add(sgens,x);
    od;
    return(Reversed(sgens));
end;

IsRigidSet := function (lst);
	if (not ForAll(lst, IsInt)) then
		return(false);
	fi;
	if (not IsSet(Reversed(lst))) then
		return(false);
	fi;
	return(true);
end;

rigid_lnc:= function(x,dim)
	# return the rigid commutators corresponding to the set x 
	#in tems of the generators of the Sylow 2-subgroup of Sym(2^)
    local n, s, lst;
	if not IsRigidSet(x) then return(fail); fi;
	if x=[] then return( () ); fi;
    s := Sigma_gens(dim);
    lst := s{x};
    return(LeftNormedComm(lst));
end;



CommutatorOfTwoRigidSets := function (r1,r2)
	# return the commutator of the two arguments as a rigid set.
	local prefix;
	if not (IsRigidSet(r1) and IsRigidSet(r2)) then
		Error("arguments must be two rigid commutators in CommutatorOfTwoRigidCommutators\n");
		return(fail);
	fi;
	if ((r1=[]) or (r2=[])) then 
		return ([]); 
	fi;
	if r1[1] = r2[1] then
		return([]);
	fi;
	if r1[1] < r2[1] then
		return (CommutatorOfTwoRigidSets (r2,r1));
	fi;
	if r2[1] in r1 then
		return([]);
	fi;
	prefix := Filtered (r1 , x -> (x > r2[1]));
	return ( Reversed(Union (prefix, [r2[1]], Intersection(r1,r2) ) ) );
end;

LeftNormedCommOfRigidSets := function(lst)
	# return the left normed commutator of the arguments in the list as a rigid set.
	# lst has to be a list of rigid sets 
	local tmplst;
	if (lst = []) then return([]); fi;
	if ( (not IsList(lst)) or (not ForAll (lst , IsRigidSet)) ) then
		Error ("arguments must be a possibly empty list of rigid sets\n");
	fi;
	if Length(lst)=1 then return (lst[1]); fi;
	tmplst:= [ CommutatorOfTwoRigidSets ( lst[1],lst[2] )];
	if Length(lst)=2 then return ( tmplst[1] ); fi;
	Append(tmplst, lst{[3..Length(lst)]});
	return(LeftNormedCommOfRigidSets(tmplst)); 
end;

dim:=8;
# we are working in Sym(2^dim)

rcombs:=ListX(Combinations([1..dim]), Reversed);;
# all rigid subsets of [1..dim]


t:=[]; 
# will contain the generators, as rigid sets, of a regular elementary abelian subgroup of Sym(2^dim)

for i in [1..dim] do
	Add(t , Reversed([1..i]));
od;

sigma:=Group(Sigma_gens(dim));
# the Sylow 2-sugroup of Sym(2^n)


rcombs_part:=Difference(rcombs,t);;
rcombs_part:=Difference(rcombs_part,[[]]);;

# Construction of the rigid generators for each of the term of normalizer chain

ngens:=[Union([[[]],Set(t)])];;

i:=-1;
while rcombs_part <> [] do
	i:=i+1;
	lst := Union (ngens);;
	new:=[];
	for x in rcombs_part do 
		if ForAll ( lst , gen -> (CommutatorOfTwoRigidSets (x , gen ) in lst) ) then
			Add(new,x);;
		fi;
	od;
	Print("The size of the ",i, "-th term of the normalizer chain is 2^", 2^dim-Size(rcombs_part)-1+Size(new),"\n");
	Print("This term is obtained by adding the following rigid commutators to the generators of the previous one:\n");
	Print(new,"\n\n");
	Add(ngens,new);;
	rcombs_part:=Difference(rcombs_part,new);;
od;


lst := [];
i:=-2;
diff:=0;
old:=dim;

Print("\n","The following table contains respectively in each column:", "\n\t" ,"1) The number i of the normalizer N_i,","\n\t",  "2) The list of the dimensions of the intersections of N_i with S_dim, ... , S_1", "\n\t", "3) log_2(Size(N_i))", "\n\t", "4) log_2(Size(N_i)/Size(N_(i-1)))", "\n\n" );

for piece in ngens do
	i:=i+1;
	piece:=Difference(piece ,[[]]);;
	lst:=Union(lst,piece);
	sizes := Reversed(List(Collected(List(lst , x -> x[1] )) , y->y[2]));;
	if (i>0) then
		old:=sum;
	fi;
	sum:=Sum(sizes);;
	if (i >=0) then Print (i,"-th term \t" ,sizes," \t-\t ", sum ," \t-\t ", sum - old ,"\n"); fi;
od;

		