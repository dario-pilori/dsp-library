clear
close all
clc

% Script to run all tests
import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.CodeCoveragePlugin;

suite = TestSuite.fromPackage('tests');
runner = TestRunner.withTextOutput;
runner.addPlugin(CodeCoveragePlugin.forFolder(pwd))
result = runner.run(suite);
