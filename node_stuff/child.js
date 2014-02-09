
var spawn = require('child_process').spawn;
 
var child = spawn('tail', ['-f', 'file.txt']);
child.stdout.on('data', function(data) {
    console.log('stdout: ' + data);
    child.kill();
});