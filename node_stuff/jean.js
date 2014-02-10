
var exec = require('child_process').execFile;
 
// var child = spawn('tail', ['-f', 'file.txt']);
// child.stdout.on('data', function(data) {
//     console.log('stdout: ' + data);
//     child.kill();
// });

var fun =function(){
   console.log("jeans running");
   //./test_exe --help
   exec('./jeans', function(err, data) {  
        console.log(err)
        console.log(data.toString());                       
    });  
}
fun();