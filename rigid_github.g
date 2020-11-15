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

Sigma_gens_as_rigid_sets := function(n)
	# Return a set of generators for Sylow 2-subgroup of Sym(2^) in the form of rigid sets
	return( List([1..n] , x ->[x]) );
end;


IsRigidSet := function (lst);
	# Check if a list of integers has strictly decreasing entries
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

NormalClosureRigid := function (bigrigidgens,smallrigidgens)
	# compute the set of rigid set generators of the normal closure
	# of the group generated by the rigid commutators corresponding to 
	# the set of rigid
	local rcombs, t, c1, c2, c3, newrigid, oldrigid;
	oldrigid:=ShallowCopy(bigrigidgens);
	newrigid:=ShallowCopy(smallrigidgens);
	for c1 in newrigid do
		for c2 in oldrigid do
			c3:=CommutatorOfTwoRigidSets(c1,c2);
			if not (c3 in newrigid) then
				Add(newrigid, c3);
			fi;
		od;
	od;
	newrigid := Difference(Set(newrigid),[[]]);;
	return (newrigid);
end;


NormalClosureChainRigid := function (dim)
	# Return the fastest subnormal chain in Sigma_n Starting from a regular abelian subgroup
	local t, c1, c2, c3, newrigid, oldrigid, sizes, i;
	t:=[]; 
	for i in [1..dim] do
		Add(t , Reversed([1..i]));
	od;
	oldrigid:=Sigma_gens_as_rigid_sets(dim);
	oldrigid:=ShallowCopy(NormalClosureRigid(oldrigid,oldrigid));
	newrigid:=ShallowCopy(t);
	sizes:=[];
	while oldrigid <> newrigid do
		Add(sizes, Length(oldrigid));
		oldrigid := NormalClosureRigid(oldrigid, newrigid);
	od;
	Add(sizes, Length(oldrigid));
	return (sizes);
end;

NormalizerRigid := function (set1, set2)
	# Returns the set of rigid generators of the normalizer of <set2> in <set1> 
	# where set1 and set 2 are saturated set of rigid sets
	local out, x, y, extset2;
	out:=[];
	extset2 := Union(set2, [[]]);
	for x in set1 do
		if ForAll (set2, y->(CommutatorOfTwoRigidSets(x,y) in extset2 ) ) then
			Add (out, x);
		fi;
	od;
	return Set(out);
end;

NormalizerRigidNonTrivial := function (set1, set2)
	# Returns the set of nontrivial rigid generators of the normalizer of <set2> in <set1> 
	# where set1 and set 2 are saturated set of rigid sets
	local out, x, y, extset2;
	out:=[];
	extset2 := Union(set2, [[]]);
	for x in set1 do
		if x <> [] then 
			if ForAll (set2, y->(CommutatorOfTwoRigidSets(x,y) in extset2 ) ) then
				Add (out, x);
			fi;
		fi;
	od;
	return Set(out);
end;

		
NormalizerRigidDifference := function (set1, set2)
	# Returns the difference with set1 of set of generators of the normalizer of <set2> in <set1> 
	# where set1 and set 2 are saturated set of rigid sets 
	local out, x, y, extset2, diff;
	out:=[];
	extset2 := Union(set2, [[]]);
	diff:=Difference(set1,set2);
	for x in diff do
		if ForAll (set2, y->(CommutatorOfTwoRigidSets(x,y) in extset2 ) ) then
			Add (out, x);
		fi;
	od;
	return Set(out);
end;

SaturationOfSetOfRigidCommutators := function (set)
	## returns the smallest saturated set of rigid commutators containing the set "set" of rigid commutators
	local inset, outset;
	inset := Set(ShallowCopy(set));
	outset:=Union( inset, Set( ListX (inset, inset, CommutatorOfTwoRigidSets) ));
	while outset <> inset do
		inset := outset;
		outset:=Union( inset, Set( ListX (inset, inset, CommutatorOfTwoRigidSets) ));
	od;
	return(outset);
end;

RigidDecomposition := function(g,dim,args...)
	# decompose a permutation as a product or rigid commutators
	local out,start,g1,g2,lp1,lp2,tmplst,x;
	if g=() then
		return [[]];
	fi;
	if Length(args) > 0 then
		start := args[1];
	else
		start:=1;
	fi;
	out:=[];
	if (1^g > 2^(dim-1)) then
		g:=Sigma_gens(dim)[1]*g;
		out:=[[start]];
#		Print(start," ",g,"\n");
	fi;
	start := start + 1;
	lp1:=ListPerm(g,2^dim);
	lp2:=ListPerm(g^Sigma_gens(dim)[1],2^dim);
	g1:=PermList(lp1{[1..2^(dim-1)]});
	g2:=PermList(lp2{[1..2^(dim-1)]});
	g1:=g1*g2;
	if g1 <> () then Append(out, RigidDecomposition(g1,dim-1,start)); fi;
	if g2 <> () then
		tmplst:=RigidDecomposition(g2,dim-1,start);
		for x in tmplst do
			Add(x,start-1);
		od;
		Append(out, tmplst ); 
	fi;
	return out;
end;

RigidCollection := function(rigidlist)
	# return an ordered list of rigid commutators having the same product as the 
	# elements in the list of rigid commutators rigidlist
	local flag, tmplst, tmprigidlst, tmpelm, pos,com;
	tmprigidlst:=DifferenceLists(rigidlist, [[]]);
	flag:=true;
	while flag do
		tmplst:= [1..(Length(tmprigidlst)-1)];
		pos:=PositionProperty(tmplst, x->(tmprigidlst[x+1] <= tmprigidlst[x]));
		if pos<>fail and tmprigidlst[pos+1]<>tmprigidlst[pos] then
			com:=CommutatorOfTwoRigidSets(tmprigidlst[pos],tmprigidlst[pos+1]);
			tmpelm := tmprigidlst[pos];
			tmprigidlst[pos] :=tmprigidlst[pos+1];
			tmprigidlst[pos+1]:= tmpelm;
			Add(tmprigidlst,com, pos+2);
		elif pos<>fail and tmprigidlst[pos+1]=tmprigidlst[pos] then
			Remove(tmprigidlst,pos+1);
			Remove(tmprigidlst,pos);
		else
			flag := false;
		fi;
	od;
	return tmprigidlst;
end;
	


dim:=11;
# we are working in Sym(2^dim) change dim accordingly to your choice

##outfile:=Concatenation("rigid_out_", String(dim), ".txt");
str := "";; outfile := OutputTextString(str,true);;


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
PrintTo(outfile,"Normalizer chain and subnormal chain of T in the Sylow 2-Subgroup of Sym(2^",dim, ")","\n\n");
while rcombs_part <> [] do
	i:=i+1;
	lst := Union (ngens);;
	new:=[];
	for x in rcombs_part do 
		if ForAll ( lst , gen -> (CommutatorOfTwoRigidSets (x , gen ) in lst) ) then
			Add(new,x);;
		fi;
	od;
	AppendTo(outfile ,"The size of the ",i, "-th term of the normalizer chain is 2^", 2^dim-Size(rcombs_part)-1+Size(new),"\n");
	AppendTo(outfile ,"This term is obtained by adding the following rigid commutators to the generators of the previous one:\n");
	AppendTo(outfile ,new,"\n\n");
	Add(ngens,new);;
	rcombs_part:=Difference(rcombs_part,new);;
od;


lst := [];
i:=-2;
diff:=0;
old:=dim;

AppendTo(outfile ,"\n","The following table contains respectively in each column:", "\n\t" ,"1) The number i of the normalizer N_i,","\n\t",  "2) The list of the dimensions of the intersections of N_i with S_",dim,", ... , S_1", "\n\t", "3) log_2(Size(N_i))", "\n\t", "4) log_2(Size(N_i)/Size(N_(i-1)))", "\n\n" );

for piece in ngens do
	i:=i+1;
	piece:=Difference(piece ,[[]]);;
	lst:=Union(lst,piece);
	sizes := Reversed(List(Collected(List(lst , x -> x[1] )) , y->y[2]));;
	if (i>0) then
		old:=sum;
	fi;
	sum:=Sum(sizes);;
	if (i >=0) then AppendTo(outfile ,i,"-th term \t" ,sizes," \t-\t ", sum ," \t-\t ", sum - old ,"\n"); fi;
od;

Print(str);

# lst:=NormalClosureChainRigid(dim);
# AppendTo(outfile ,"\n", "The sunbnormality defect of t is ",Length(lst)," and subgroup of the chain have orders");
# for i in lst do
	# AppendTo(outfile ," 2^", i);
# od;
# AppendTo(outfile ,"\n");

# tgroup:=Group(List(t,x->rigid_lnc(x,dim)));
# u_rigid_gens:=NormalizerRigidNonTrivial(rcombs,t);
# ugroup:=Group(List(u_rigid_gens,x->rigid_lnc(x,dim)));
# n1_rigid_gens:=NormalizerRigidNonTrivial(rcombs,u_rigid_gens);
# n1group:= Group(List(n1_rigid_gens,x->rigid_lnc(x,dim)));


# sym:=SymmetricGroup(2^dim);
# n1:=Normalizer(sym,ugroup);
# time;

# n1=n1group;

