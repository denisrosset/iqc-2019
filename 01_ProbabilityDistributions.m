
%% Exported from Jupyter Notebook
% Run each section by placing your cursor in it and pressing Ctrl+Enter

%% Markdown Cell:
% Before we start, let us ensure that RepLAB is on the path. Please modify the line below to correspond to your RepLAB installation folder (or add RepLAB to your MATLAB/Octave path).

%% Code Cell[11]:

run ~/software/replab/replab_addpaths(2,1)
install_sdpt3

%% Markdown Cell:
% # Permutation groups, application to probability distributions
% 
% Reading recommendation: [Linear Representations of Finite Groups](https://link.springer.com/book/10.1007/978-1-4684-9458-7), Jean-Pierre Serre
% 
% Consider the probability distribution $P_\text{A}(a)$ for $a=1,2,3,4$ corresponding to the observations of a measurement with four outcomes.
% 
% In the *device-independent* framework [see the review](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.86.419), the particular labels attached to these outcomes do not matter.
% 
% Thus, the index $a$ can be relabeled without changing the underlying physics using a [permutation](https://en.wikipedia.org/wiki/Permutation).
% 
% For permutations, we use the convention that permutations are described using a row vector of their images; indices are 1-based (Matlab/Octave).
% 
% Thus we write as follows 
% 
% - a cyclic shift of all 4 elements,
% - a cyclic shift of the first three elements.

%% Code Cell[12]:

c4 = [2 3 4 1]
c3 = [2 3 1 4]

%% Markdown Cell:
% Remark that when writing permutations, we implicitly assume the size of the domain (here, 4). In RepLAB, permutations of different sizes cannot be composed. This is in contrast to most computer algebra systems that represent permutations in a more abstract fashion, using a product of cycles.
% 
% To write the composition $c_4 \cdot c_3$, we need to decide whether the action of $c_4$ and $c_3$ on labels is a [left action or a right action](https://en.wikipedia.org/wiki/Group_action_(mathematics)#Definition).
% 
% In physics, the left action convention is prevalent, while in computational group theory, the right action is used.
% 
% RepLAB uses the left action convention.
% 
% To compose two permutations, we can construct the group of all permutations of four elements, $\mathcal{S}_4$, and use the `compose` method of that group.
% 
% We then verify that the left action convention works.

%% Code Cell[14]:

S4 = replab.S(4) % shortcut for S4 = replab.Permutations(4)
c4c3 = S4.compose(c4, c3) % it corresponds to c4(c3)
c4inv = S4.inverse(c4)

%% Code Cell[15]:

i = 1 % replace by 2,3,4 as well
img1 = c4c3(i)
img2 = c4(c3(i))

%% Markdown Cell:
% ## Representations
% 
% Linear representations represent the action of a group on a vector space, in a way that preserves the group structure and the vector space structure, see [Definition](https://en.wikipedia.org/wiki/Representation_theory_of_finite_groups#Linear_representations).
% 
% A permutation group $G$ that permutes $n$ indices (i.e. $G \subseteq \mathcal{S}_4$) acts naturally on vectors $\in \mathbb{R}^n$ by permuting their coordinates; this is the *natural representation*.
% 
% Finite groups are fully described by their [set of generators](https://en.wikipedia.org/wiki/Generating_set_of_a_group), with the property that every group element can be written as a product of generators.
% 
% Thus, we can describe a representation by the matrices that represent the action of the generators.
% 
% Compare below `S4.generator(1)`, `S4.generator(2)` to `rho.images{1}` and `rho.images{2}` respectively.

%% Code Cell[16]:

S4
rho = S4.naturalRep

%% Markdown Cell:
% ### Representation: image of a group element
% 
% A (unitary) representation provides a unitary matrix $\rho_g$ for each group element $g$.
% 
% We can ask RepLAB to compute the image of any group element.

%% Code Cell[18]:

g = [2 4 1 3];
full(rho.image(g))% "full" used to show the full matrix instead of the sparse matrix elements

%% Markdown Cell:
% ### Decomposing representations
% 
% Are there vector subspaces of $\mathcal{R}^4$ that are invariant under $\rho_g$ for any $g\in G$?
% 
% For sure, the vector $\vec{v} = (1,1,1,1)^\top$ is invariant under any relabeling -- and thus the subspace spanned by $\vec{v}$ is invariant as well.
% 
% We can thus restrict $\rho$ to that subspace; the method `subRep` requires row vectors. Per our convention, those row vectors need to be orthogonal but they do not need to be normalized.
% 
% We then look at how this subrepresentation acts on its subspace.

%% Code Cell[19]:

rho1 = rho.subRep([1 1 1 1])
rho1.image(S4.generator(1))
rho1.image(S4.generator(2))

%% Markdown Cell:
% The image of any group element is the identity for $\rho_1$: by convention, such representations are called *trivial representations*.
% 
% The subspace orthogonal to $\vec{v}$ is also a subrepresentation (see [Maschke's theorem](https://en.wikipedia.org/wiki/Maschke%27s_theorem), and Jean-Pierre Serre chapter 1). For the symmetric group, that subspace corresponds to the [standard representation](https://en.wikipedia.org/wiki/Representation_theory_of_the_symmetric_group#Low-dimensional_representations).

%% Code Cell[25]:

rho2 = rho.subRep([1 -1 1 -1; 1 1 -1 -1; 1 -1 -1 1])

%% Code Cell[26]:

replab.RepLaws(rho2).check

%% Markdown Cell:
% ### Decomposing representations
% 
% To decompose a representation, one finds invariant subspaces and splits subrepresentations accordingly.
% 
% RepLAB provides a way to compute this decomposition numerically.
% 
% We do that, and verify that the subrepresentations match our intuition above.

%% Code Cell[27]:

I = rho.decomposition

%% Code Cell[29]:

I.component(1).copy(1)

%% Code Cell[30]:

I.component(2).copy(1)

%% Code Cell[43]:

[C Crep] = replab.quantum.clifford_qudit(3)

%% Code Cell[47]:

UxU = kron(Crep, conj(Crep))

%% Code Cell[48]:

I=UxU.decomposition

%% Code Cell[49]:

I.component(1).copy(1)

%% Code Cell[ ]:


