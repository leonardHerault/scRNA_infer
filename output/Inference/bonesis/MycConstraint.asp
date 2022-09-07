#!/usr/bin/env bash
clingo -W no-atom-undefined -t 16 0 --single-shot --project -c bounded_nonreach=0 "${@}" - <<EOF
#program base.
{clause(N,1..C,L,S): in(L,N,S), maxC(N,C), node(N), node(L)}.
:- clause(N,_,L,S), clause(N,_,L,-S).
1 { constant(N,(-1;1)) } 1 :- node(N), not clause(N,_,_,_).
constant(N) :- constant(N,_).
size(N,C,X) :- X = #count {L,S: clause(N,C,L,S)}; clause(N,C,_,_); maxC(N,_).
:- clause(N,C,_,_); not clause(N,C-1,_,_); C > 1; maxC(N,_).
:- size(N,C1,X1); size(N,C2,X2); X1 < X2; C1 > C2; maxC(N,_).
:- size(N,C1,X); size(N,C2,X); C1 > C2; mindiff(N,C1,C2,L1) ; mindiff(N,C2,C1,L2) ; L1 < L2; maxC(N,_).
clausediff(N,C1,C2,L) :- clause(N,C1,L,_);not clause(N,C2,L,_);clause(N,C2,_,_), C1 != C2; maxC(N,_).
mindiff(N,C1,C2,L) :- clausediff(N,C1,C2,L); L <= L' : clausediff(N,C1,C2,L'), clause(N,C1,L',_), C1!=C2; maxC(N,_).
:- size(N,C1,X1); size(N,C2,X2); C1 != C2; X1 <= X2; clause(N,C2,L,S) : clause(N,C1,L,S); maxC(N,_).
nbnode(15).
node("Bclaf1").
node("CDK46CycD").
node("Myc").
node("Cebpa").
node("Gata2").
node("Egr1").
node("Junb").
node("CIPKIP").
node("Fli1").
node("Zfpm1").
node("Gata1").
node("Klf1").
node("Ikzf1").
node("Tal1").
node("Spi1").
in("Bclaf1","CDK46CycD",1).
in("Bclaf1","Bclaf1",1).
in("Bclaf1","Myc",1).
in("Bclaf1","Cebpa",1).
in("Myc","Junb",1).
in("Myc","CDK46CycD",1).
in("Myc","Bclaf1",1).
in("Myc","Myc",1).
in("Myc","Cebpa",1).
in("Cebpa","Gata2",1).
in("Cebpa","Myc",1).
in("Cebpa","Cebpa",1).
in("Cebpa","Spi1",1).
in("Gata2","Ikzf1",1).
in("Gata2","Zfpm1",1).
in("Gata2","Fli1",1).
in("Gata2","Cebpa",1).
in("Gata2","Egr1",1).
in("Gata2","Tal1",1).
in("Gata2","Gata2",1).
in("Gata2","Gata1",1).
in("Gata2","Spi1",-1).
in("Egr1","Junb",1).
in("Egr1","Egr1",1).
in("Egr1","Gata2",1).
in("Egr1","CIPKIP",1).
in("Junb","Gata2",1).
in("Junb","CDK46CycD",1).
in("Junb","Egr1",1).
in("Junb","Junb",1).
in("Junb","Fli1",1).
in("Junb","Myc",1).
in("Junb","CIPKIP",1).
in("Fli1","Fli1",1).
in("Fli1","Zfpm1",1).
in("Fli1","Gata1",1).
in("Fli1","Klf1",-1).
in("Zfpm1","Cebpa",-1).
in("Zfpm1","Gata2",-1).
in("Gata1","Klf1",1).
in("Gata1","Gata1",1).
in("Gata1","Cebpa",-1).
in("Gata1","Fli1",1).
in("Gata1","Gata2",-1).
in("Gata1","Ikzf1",-1).
in("Gata1","Spi1",-1).
in("Gata1","Tal1",1).
in("Gata1","Zfpm1",1).
in("Klf1","Gata1",1).
in("Klf1","Fli1",-1).
in("Ikzf1","Ikzf1",1).
in("Ikzf1","Cebpa",-1).
in("Ikzf1","Gata1",-1).
in("Tal1","Cebpa",-1).
in("Spi1","Cebpa",1).
in("Spi1","Myc",1).
in("Spi1","Gata1",-1).
in("Spi1","Gata2",-1).
in("Spi1","Spi1",1).
in("Spi1","Tal1",-1).
maxC("Bclaf1",2).
maxC("CDK46CycD",3).
maxC("Myc",3).
maxC("Cebpa",3).
maxC("Gata2",3).
maxC("Egr1",3).
maxC("Junb",3).
maxC("CIPKIP",2).
maxC("Fli1",3).
maxC("Zfpm1",3).
maxC("Gata1",3).
maxC("Klf1",2).
maxC("Ikzf1",3).
maxC("Tal1",3).
maxC("Spi1",3).
obs("iHSC","Egr1",-1).
obs("iHSC","Junb",-1).
obs("iHSC","Bclaf1",1).
obs("iHSC","Myc",-1).
obs("iHSC","Fli1",1).
obs("iHSC","Gata2",1).
obs("iHSC","Spi1",-1).
obs("iHSC","Cebpa",-1).
obs("iHSC","Gata1",-1).
obs("iHSC","Klf1",-1).
obs("iHSC","Tal1",1).
obs("iHSC","Ikzf1",-1).
obs("iHSC","CDK46CycD",-1).
obs("iHSC","CIPKIP",-1).
obs("srHSC","Junb",-1).
obs("srHSC","Bclaf1",1).
obs("srHSC","Myc",-1).
obs("srHSC","Fli1",1).
obs("srHSC","Gata2",-1).
obs("srHSC","Spi1",-1).
obs("srHSC","Cebpa",-1).
obs("srHSC","Gata1",-1).
obs("srHSC","Klf1",-1).
obs("srHSC","CDK46CycD",1).
obs("srHSC","CIPKIP",-1).
obs("qHSC","Egr1",1).
obs("qHSC","Junb",1).
obs("qHSC","Bclaf1",-1).
obs("qHSC","Myc",1).
obs("qHSC","Fli1",1).
obs("qHSC","Gata2",1).
obs("qHSC","Spi1",-1).
obs("qHSC","Cebpa",-1).
obs("qHSC","Gata1",-1).
obs("qHSC","Klf1",-1).
obs("qHSC","Tal1",1).
obs("qHSC","Ikzf1",-1).
obs("qHSC","CDK46CycD",1).
obs("qHSC","CIPKIP",1).
obs("preDiff","Egr1",-1).
obs("preDiff","Junb",-1).
obs("preDiff","Bclaf1",1).
obs("preDiff","Myc",1).
obs("preDiff","Fli1",-1).
obs("preDiff","Gata2",-1).
obs("preDiff","Spi1",1).
obs("preDiff","Cebpa",-1).
obs("preDiff","Gata1",-1).
obs("preDiff","Klf1",-1).
obs("preDiff","CDK46CycD",-1).
obs("preDiff","CIPKIP",-1).
obs("pLymph","Egr1",-1).
obs("pLymph","Junb",-1).
obs("pLymph","Myc",-1).
obs("pLymph","Fli1",-1).
obs("pLymph","Gata2",1).
obs("pLymph","Spi1",1).
obs("pLymph","Cebpa",-1).
obs("pLymph","Gata1",-1).
obs("pLymph","Klf1",-1).
obs("pLymph","Tal1",-1).
obs("pLymph","Ikzf1",1).
obs("pLymph","CDK46CycD",-1).
obs("pNeuMast","Egr1",-1).
obs("pNeuMast","Junb",-1).
obs("pNeuMast","Bclaf1",-1).
obs("pNeuMast","Fli1",-1).
obs("pNeuMast","Spi1",1).
obs("pNeuMast","Cebpa",1).
obs("pNeuMast","Gata1",-1).
obs("pNeuMast","Klf1",-1).
obs("pNeuMast","Tal1",-1).
obs("pNeuMast","Zfpm1",-1).
obs("pNeuMast","CDK46CycD",-1).
obs("pNeuMast","CIPKIP",-1).
obs("pMk","Egr1",-1).
obs("pMk","Junb",-1).
obs("pMk","Bclaf1",-1).
obs("pMk","Fli1",1).
obs("pMk","Gata2",-1).
obs("pMk","Spi1",-1).
obs("pMk","Cebpa",-1).
obs("pMk","Gata1",1).
obs("pMk","Klf1",-1).
obs("pMk","Tal1",1).
obs("pMk","Ikzf1",-1).
obs("pMk","Zfpm1",1).
obs("pMk","CIPKIP",-1).
obs("pEr","Egr1",-1).
obs("pEr","Junb",-1).
obs("pEr","Bclaf1",-1).
obs("pEr","Fli1",-1).
obs("pEr","Gata2",-1).
obs("pEr","Spi1",-1).
obs("pEr","Cebpa",-1).
obs("pEr","Gata1",1).
obs("pEr","Klf1",1).
obs("pEr","Zfpm1",1).
obs("pEr","CDK46CycD",-1).
obs("pEr","CIPKIP",-1).
obs("G0pMk","Bclaf1",-1).
obs("G0pMk","Myc",-1).
obs("G0pMk","Fli1",1).
obs("G0pMk","Gata2",-1).
obs("G0pMk","Spi1",-1).
obs("G0pMk","Cebpa",-1).
obs("G0pMk","Gata1",1).
obs("G0pMk","Klf1",-1).
obs("G0pMk","Tal1",1).
obs("G0pMk","Ikzf1",-1).
obs("G0pMk","Zfpm1",1).
obs("G0pMk","CIPKIP",1).
obs("G2MpNeuMast","Egr1",-1).
obs("G2MpNeuMast","Junb",-1).
obs("G2MpNeuMast","Fli1",-1).
obs("G2MpNeuMast","Spi1",1).
obs("G2MpNeuMast","Cebpa",1).
obs("G2MpNeuMast","Gata1",-1).
obs("G2MpNeuMast","Klf1",-1).
obs("G2MpNeuMast","Tal1",-1).
obs("G2MpNeuMast","Ikzf1",-1).
obs("G2MpNeuMast","Zfpm1",-1).
obs("G2MpNeuMast","CDK46CycD",1).
obs("G2MpNeuMast","CIPKIP",-1).
obs("zeros","Egr1",-1).
obs("zeros","Junb",-1).
obs("zeros","Tal1",-1).
obs("zeros","Bclaf1",-1).
obs("zeros","Myc",-1).
obs("zeros","Fli1",-1).
obs("zeros","Gata2",-1).
obs("zeros","Ikzf1",-1).
obs("zeros","Spi1",-1).
obs("zeros","Cebpa",-1).
obs("zeros","Gata1",-1).
obs("zeros","Klf1",-1).
obs("zeros","Zfpm1",-1).
obs("zeros","CIPKIP",-1).
obs("zeros","CDK46CycD",-1).
obs("nonPrimed","Ikzf1",-1).
obs("nonPrimed","Spi1",-1).
obs("nonPrimed","Cebpa",-1).
obs("nonPrimed","Gata1",-1).
obs("nonPrimed","Klf1",-1).
1 {cfg(X,N,(-1;1))} 1 :- cfg(X), node(N).
cfg(X,N,V) :- bind_cfg(X,O), obs(O,N,V), node(N).
eval(X,N,C,-1) :- clause(N,C,L,-V), mcfg(X,L,V), not clamped(X,N,_).
eval(X,N,C,1) :- mcfg(X,L,V): clause(N,C,L,V); clause(N,C,_,_), mcfg(X,_,_), not clamped(X,N,_).
eval(X,N,1) :- eval(X,N,C,1), clause(N,C,_,_).
eval(X,N,-1) :- eval(X,N,C,-1): clause(N,C,_,_); clause(N,_,_,_), mcfg(X,_,_).
eval(X,N,V) :- clamped(X,N,V).
eval(X,N,V) :- constant(N,V), mcfg(X,_,_), not clamped(X,N,_).
mcfg(X,N,V) :- ext(X,N,V).
cfg(__bocfg22,N,-1); cfg(__bocfg22,N,1) :- node(N).
cfg(__bocfg22,N,-V) :- cfg(__bocfg22,N,V), saturate(__bocfg22).
saturate(__bocfg22) :- valid(__bocfg22,Z): expect_valid(__bocfg22,Z).
:- not saturate(__bocfg22).
expect_valid(__bocfg22,__bocond23).
cfg(X,N,V) :- bind_cfg(X,O,mutant(M)), obs(O,N,V), node(N), not mutant(M,N,_).
cfg(X,N,V) :- bind_cfg(X,O,mutant(M)), obs(O,_,_), node(N), mutant(M,N,V), not weak_mutant(M,N,V).
cfg(X,N,V) :- bind_cfg(X,O,mutant(M)), obs(O,N,V), node(N), mutant(M,N,W), weak_mutant(M,N,W).
expect_valid(__bocfg22,__bocond26).
expect_valid(__bocfg22,__bocond29).
obs("pLymph").
cfg("pLymph").
bind_cfg("pLymph","pLymph").
mcfg(__bofp0,N,V) :- cfg("pLymph",N,V).
:- cfg("pLymph",N,V), eval(__bofp0,N,-V).
obs("pEr").
cfg("pEr").
bind_cfg("pEr","pEr").
mcfg(__bofp1,N,V) :- cfg("pEr",N,V).
:- cfg("pEr",N,V), eval(__bofp1,N,-V).
obs("pMk").
cfg("pMk").
bind_cfg("pMk","pMk").
mcfg(__bofp2,N,V) :- cfg("pMk",N,V).
:- cfg("pMk",N,V), eval(__bofp2,N,-V).
obs("pNeuMast").
cfg("pNeuMast").
bind_cfg("pNeuMast","pNeuMast").
mcfg(__bofp3,N,V) :- cfg("pNeuMast",N,V).
:- cfg("pNeuMast",N,V), eval(__bofp3,N,-V).
obs("iHSC").
cfg("iHSC").
bind_cfg("iHSC","iHSC").
obs("srHSC").
cfg("srHSC").
bind_cfg("srHSC","srHSC").
mcfg(__boreach4,N,V) :- cfg("iHSC",N,V).
ext(__boreach4,N,V) :- eval(__boreach4,N,V), cfg("srHSC",N,V).
{ext(__boreach4,N,V)} :- eval(__boreach4,N,V), cfg("srHSC",N,-V).
:- cfg("srHSC",N,V), not mcfg(__boreach4,N,V).
:- cfg("srHSC",N,V), ext(__boreach4,N,-V), not ext(__boreach4,N,V).
obs("qHSC").
cfg("qHSC").
bind_cfg("qHSC","qHSC").
mcfg(__boreach5,N,V) :- cfg("iHSC",N,V).
ext(__boreach5,N,V) :- eval(__boreach5,N,V), cfg("qHSC",N,V).
{ext(__boreach5,N,V)} :- eval(__boreach5,N,V), cfg("qHSC",N,-V).
:- cfg("qHSC",N,V), not mcfg(__boreach5,N,V).
:- cfg("qHSC",N,V), ext(__boreach5,N,-V), not ext(__boreach5,N,V).
mcfg(__boreach6,N,V) :- cfg("iHSC",N,V).
ext(__boreach6,N,V) :- eval(__boreach6,N,V), cfg("pLymph",N,V).
{ext(__boreach6,N,V)} :- eval(__boreach6,N,V), cfg("pLymph",N,-V).
:- cfg("pLymph",N,V), not mcfg(__boreach6,N,V).
:- cfg("pLymph",N,V), ext(__boreach6,N,-V), not ext(__boreach6,N,V).
obs("preDiff").
cfg("preDiff").
bind_cfg("preDiff","preDiff").
mcfg(__boreach7,N,V) :- cfg("iHSC",N,V).
ext(__boreach7,N,V) :- eval(__boreach7,N,V), cfg("preDiff",N,V).
{ext(__boreach7,N,V)} :- eval(__boreach7,N,V), cfg("preDiff",N,-V).
:- cfg("preDiff",N,V), not mcfg(__boreach7,N,V).
:- cfg("preDiff",N,V), ext(__boreach7,N,-V), not ext(__boreach7,N,V).
mcfg(__boreach8,N,V) :- cfg("preDiff",N,V).
ext(__boreach8,N,V) :- eval(__boreach8,N,V), cfg("pEr",N,V).
{ext(__boreach8,N,V)} :- eval(__boreach8,N,V), cfg("pEr",N,-V).
:- cfg("pEr",N,V), not mcfg(__boreach8,N,V).
:- cfg("pEr",N,V), ext(__boreach8,N,-V), not ext(__boreach8,N,V).
mcfg(__boreach9,N,V) :- cfg("iHSC",N,V).
ext(__boreach9,N,V) :- eval(__boreach9,N,V), cfg("preDiff",N,V).
{ext(__boreach9,N,V)} :- eval(__boreach9,N,V), cfg("preDiff",N,-V).
:- cfg("preDiff",N,V), not mcfg(__boreach9,N,V).
:- cfg("preDiff",N,V), ext(__boreach9,N,-V), not ext(__boreach9,N,V).
mcfg(__boreach10,N,V) :- cfg("preDiff",N,V).
ext(__boreach10,N,V) :- eval(__boreach10,N,V), cfg("pMk",N,V).
{ext(__boreach10,N,V)} :- eval(__boreach10,N,V), cfg("pMk",N,-V).
:- cfg("pMk",N,V), not mcfg(__boreach10,N,V).
:- cfg("pMk",N,V), ext(__boreach10,N,-V), not ext(__boreach10,N,V).
mcfg(__boreach11,N,V) :- cfg("iHSC",N,V).
ext(__boreach11,N,V) :- eval(__boreach11,N,V), cfg("preDiff",N,V).
{ext(__boreach11,N,V)} :- eval(__boreach11,N,V), cfg("preDiff",N,-V).
:- cfg("preDiff",N,V), not mcfg(__boreach11,N,V).
:- cfg("preDiff",N,V), ext(__boreach11,N,-V), not ext(__boreach11,N,V).
mcfg(__boreach12,N,V) :- cfg("preDiff",N,V).
ext(__boreach12,N,V) :- eval(__boreach12,N,V), cfg("pNeuMast",N,V).
{ext(__boreach12,N,V)} :- eval(__boreach12,N,V), cfg("pNeuMast",N,-V).
:- cfg("pNeuMast",N,V), not mcfg(__boreach12,N,V).
:- cfg("pNeuMast",N,V), ext(__boreach12,N,-V), not ext(__boreach12,N,V).
mcfg(__boreach13,N,V) :- cfg("srHSC",N,V).
ext(__boreach13,N,V) :- eval(__boreach13,N,V), cfg("iHSC",N,V).
{ext(__boreach13,N,V)} :- eval(__boreach13,N,V), cfg("iHSC",N,-V).
:- cfg("iHSC",N,V), not mcfg(__boreach13,N,V).
:- cfg("iHSC",N,V), ext(__boreach13,N,-V), not ext(__boreach13,N,V).
mcfg(__boreach14,N,V) :- cfg("qHSC",N,V).
ext(__boreach14,N,V) :- eval(__boreach14,N,V), cfg("iHSC",N,V).
{ext(__boreach14,N,V)} :- eval(__boreach14,N,V), cfg("iHSC",N,-V).
:- cfg("iHSC",N,V), not mcfg(__boreach14,N,V).
:- cfg("iHSC",N,V), ext(__boreach14,N,-V), not ext(__boreach14,N,V).
mcfg((__bononreach15,1..K),N,V) :- reach_steps(__bononreach15,K), cfg("preDiff",N,V).
ext((__bononreach15,I),N,V) :- eval((__bononreach15,I),N,V), not locked((__bononreach15,I),N).
reach_bad(__bononreach15,I,N) :- cfg("preDiff",N,V), cfg("qHSC",N,V), ext((__bononreach15,I),N,-V), not ext((__bononreach15,I),N,V).
locked((__bononreach15,I+1..K),N) :- reach_bad(__bononreach15,I,N), reach_steps(__bononreach15,K), I < K.
nr_ok(__bononreach15) :- reach_steps(__bononreach15,K), cfg("qHSC",N,V), not mcfg((__bononreach15,K),N,V).
:- not nr_ok(__bononreach15).
reach_steps(__bononreach15,K) :- nbnode(K), bounded_nonreach <= 0.
reach_steps(__bononreach15,bounded_nonreach) :- bounded_nonreach > 0.
mcfg((__bononreach16,1..K),N,V) :- reach_steps(__bononreach16,K), cfg("preDiff",N,V).
ext((__bononreach16,I),N,V) :- eval((__bononreach16,I),N,V), not locked((__bononreach16,I),N).
reach_bad(__bononreach16,I,N) :- cfg("preDiff",N,V), cfg("srHSC",N,V), ext((__bononreach16,I),N,-V), not ext((__bononreach16,I),N,V).
locked((__bononreach16,I+1..K),N) :- reach_bad(__bononreach16,I,N), reach_steps(__bononreach16,K), I < K.
nr_ok(__bononreach16) :- reach_steps(__bononreach16,K), cfg("srHSC",N,V), not mcfg((__bononreach16,K),N,V).
:- not nr_ok(__bononreach16).
reach_steps(__bononreach16,K) :- nbnode(K), bounded_nonreach <= 0.
reach_steps(__bononreach16,bounded_nonreach) :- bounded_nonreach > 0.
mcfg((__bononreach17,1..K),N,V) :- reach_steps(__bononreach17,K), cfg("preDiff",N,V).
ext((__bononreach17,I),N,V) :- eval((__bononreach17,I),N,V), not locked((__bononreach17,I),N).
reach_bad(__bononreach17,I,N) :- cfg("preDiff",N,V), cfg("iHSC",N,V), ext((__bononreach17,I),N,-V), not ext((__bononreach17,I),N,V).
locked((__bononreach17,I+1..K),N) :- reach_bad(__bononreach17,I,N), reach_steps(__bononreach17,K), I < K.
nr_ok(__bononreach17) :- reach_steps(__bononreach17,K), cfg("iHSC",N,V), not mcfg((__bononreach17,K),N,V).
:- not nr_ok(__bononreach17).
reach_steps(__bononreach17,K) :- nbnode(K), bounded_nonreach <= 0.
reach_steps(__bononreach17,bounded_nonreach) :- bounded_nonreach > 0.
obs("zeros").
cfg("zeros").
bind_cfg("zeros","zeros").
mcfg((__bononreach18,1..K),N,V) :- reach_steps(__bononreach18,K), cfg("zeros",N,V).
ext((__bononreach18,I),N,V) :- eval((__bononreach18,I),N,V), not locked((__bononreach18,I),N).
reach_bad(__bononreach18,I,N) :- cfg("zeros",N,V), cfg("pNeuMast",N,V), ext((__bononreach18,I),N,-V), not ext((__bononreach18,I),N,V).
locked((__bononreach18,I+1..K),N) :- reach_bad(__bononreach18,I,N), reach_steps(__bononreach18,K), I < K.
nr_ok(__bononreach18) :- reach_steps(__bononreach18,K), cfg("pNeuMast",N,V), not mcfg((__bononreach18,K),N,V).
:- not nr_ok(__bononreach18).
reach_steps(__bononreach18,1).
mcfg((__bononreach19,1..K),N,V) :- reach_steps(__bononreach19,K), cfg("zeros",N,V).
ext((__bononreach19,I),N,V) :- eval((__bononreach19,I),N,V), not locked((__bononreach19,I),N).
reach_bad(__bononreach19,I,N) :- cfg("zeros",N,V), cfg("pMk",N,V), ext((__bononreach19,I),N,-V), not ext((__bononreach19,I),N,V).
locked((__bononreach19,I+1..K),N) :- reach_bad(__bononreach19,I,N), reach_steps(__bononreach19,K), I < K.
nr_ok(__bononreach19) :- reach_steps(__bononreach19,K), cfg("pMk",N,V), not mcfg((__bononreach19,K),N,V).
:- not nr_ok(__bononreach19).
reach_steps(__bononreach19,1).
mcfg((__bononreach20,1..K),N,V) :- reach_steps(__bononreach20,K), cfg("zeros",N,V).
ext((__bononreach20,I),N,V) :- eval((__bononreach20,I),N,V), not locked((__bononreach20,I),N).
reach_bad(__bononreach20,I,N) :- cfg("zeros",N,V), cfg("pEr",N,V), ext((__bononreach20,I),N,-V), not ext((__bononreach20,I),N,V).
locked((__bononreach20,I+1..K),N) :- reach_bad(__bononreach20,I,N), reach_steps(__bononreach20,K), I < K.
nr_ok(__bononreach20) :- reach_steps(__bononreach20,K), cfg("pEr",N,V), not mcfg((__bononreach20,K),N,V).
:- not nr_ok(__bononreach20).
reach_steps(__bononreach20,1).
mcfg((__bononreach21,1..K),N,V) :- reach_steps(__bononreach21,K), cfg("zeros",N,V).
ext((__bononreach21,I),N,V) :- eval((__bononreach21,I),N,V), not locked((__bononreach21,I),N).
reach_bad(__bononreach21,I,N) :- cfg("zeros",N,V), cfg("pLymph",N,V), ext((__bononreach21,I),N,-V), not ext((__bononreach21,I),N,V).
locked((__bononreach21,I+1..K),N) :- reach_bad(__bononreach21,I,N), reach_steps(__bononreach21,K), I < K.
nr_ok(__bononreach21) :- reach_steps(__bononreach21,K), cfg("pLymph",N,V), not mcfg((__bononreach21,K),N,V).
:- not nr_ok(__bononreach21).
reach_steps(__bononreach21,1).
mcfg(__bocfg24,N,V) :- cfg(__bocfg22,N,V).
valid(__bocfg22,__bocond23) :- cfg(__bocfg22,N,V), eval(__bocfg24,N,-V).
valid(__bocfg22,__bocond23) :- cfg(__bocfg22,N,V): obs("zeros",N,V), node(N).
valid(__bocfg22,__bocond23) :- cfg(__bocfg22,N,V): obs("pEr",N,V), node(N).
valid(__bocfg22,__bocond23) :- cfg(__bocfg22,N,V): obs("pLymph",N,V), node(N).
valid(__bocfg22,__bocond23) :- cfg(__bocfg22,N,V): obs("pNeuMast",N,V), node(N).
valid(__bocfg22,__bocond23) :- cfg(__bocfg22,N,V): obs("pMk",N,V), node(N).
mcfg(__bocfg25,N,V) :- cfg("iHSC",N,V).
mcfg(__bocfg25,N,V) :- eval(__bocfg25,N,V).
valid(__bocfg22,__bocond23) :- cfg(__bocfg22,N,V), not mcfg(__bocfg25,N,V).
mutant(1,"Myc",-1).
cfg(("iHSC",0)).
bind_cfg(("iHSC",0),"iHSC",mutant(1)).
cfg(("iHSC",0),"Zfpm1",-1).
obs("nonPrimed").
mcfg(__bocfg27,N,V) :- cfg(__bocfg22,N,V).
valid(__bocfg22,__bocond26) :- cfg(__bocfg22,N,V), eval(__bocfg27,N,-V).
valid(__bocfg22,__bocond26) :- cfg(__bocfg22,N,V): obs("nonPrimed",N,V), node(N).
clamped(__bocfg27,N,V) :- mutant(1,N,V).
mcfg(__bocfg28,N,V) :- cfg(("iHSC",0),N,V).
mcfg(__bocfg28,N,V) :- eval(__bocfg28,N,V).
valid(__bocfg22,__bocond26) :- cfg(__bocfg22,N,V), not mcfg(__bocfg28,N,V).
clamped(__bocfg28,N,V) :- mutant(1,N,V).
cfg(("iHSC",1)).
bind_cfg(("iHSC",1),"iHSC",mutant(1)).
cfg(("iHSC",1),"Zfpm1",1).
mcfg(__bocfg30,N,V) :- cfg(__bocfg22,N,V).
valid(__bocfg22,__bocond29) :- cfg(__bocfg22,N,V), eval(__bocfg30,N,-V).
valid(__bocfg22,__bocond29) :- cfg(__bocfg22,N,V): obs("nonPrimed",N,V), node(N).
clamped(__bocfg30,N,V) :- mutant(1,N,V).
mcfg(__bocfg31,N,V) :- cfg(("iHSC",1),N,V).
mcfg(__bocfg31,N,V) :- eval(__bocfg31,N,V).
valid(__bocfg22,__bocond29) :- cfg(__bocfg22,N,V), not mcfg(__bocfg31,N,V).
clamped(__bocfg31,N,V) :- mutant(1,N,V).
#show clause/4.
#show constant/2.

EOF
