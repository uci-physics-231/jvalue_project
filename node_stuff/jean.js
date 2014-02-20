
// var exec = require('child_process').execFile;
 
// var child = spawn('tail', ['-f', 'file.txt']);
// child.stdout.on('data', function(data) {
//     console.log('stdout: ' + data);
//     child.kill();
// });
// function run_and_kill() {
// 	var exec = require('child_process').exec,
// 		ls = exec('./test_exe');
// 	console.log('Child process started: %d', ls.pid);
// 	
// 	ls.on('exit', function(code, signal) {
// 		console.log('exit with code %s and signal %s', code, signal);
// 	});
// 	ls.kill();
// }
// run_and_kill();
// 			
// var fun =function(){
// var exec = require('child_process').execFile;
//    console.log("jeans running");
//    //./test_exe --help
//    var shit = exec('./test_exe');
//    //,[' --help'], function(err, data) { 
// //    exec('./jeans', function(err, data) {  
//      //   console.log(err)
//      //   console.log(data.toString());                       
//    // });  
// }
// fun();

var spawn = require('child_process').spawn;
 
var child = spawn('./test_exe', ['--help']);
child.stdout.on('data', function(data) {
    console.log('stdout: ' + data);
    child.kill();
});
