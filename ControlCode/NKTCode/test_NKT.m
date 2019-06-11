%% test_NKT.m
% Script to test the basic functionality of the NKT source
% author: G. Ruane 
% created: 2019 March 22

%clear;

% addpath('/home/hcst/hcst_lib/');

bench.info.NKT_lib_PATH = '/home/hcst/hcst_lib/NKT/';

% bench.NKT.CONNECTED = false; 
% 
% disp('*** Connecting to NKT ... ***');
% 
% % Create Nkt object using nkt_mod_falco.py
% bench.NKT.nktobj = py.nkt_mod_falco.Nkt('/dev/ttyNKT',115200);
% bench.NKT.CONNECTED = true;



disp('*** NKT connected. ***');
%%
bench.NKT.nktobj.get_emission()
% bench.NKT.nktobj.set_emission(true);
%%

% tb_setUpNKT(bench);

%%

%tb_NKT_setPowerLevel(bench,powerlevel)
%tb_NKT_getPowerLevel(bench)
%tb_NKT_setEmission(bench,true)
%tb_NKT_setWvlRange(bench,573,577)
%tb_NKT_getWvlRange(bench)
%%
%tb_cleanUpNKT(bench)
