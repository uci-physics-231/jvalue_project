var http = require('http');
var url = require("url");


// var childProcess = require('child_process').spawn(),
//      ls;
//      
// // ls -l
// ls = childProcess.exec('ls -l', function (error, stdout, stderr) {
//    if (error) {
//      console.log(error.stack);
//      console.log('Error code: '+error.code);
//      console.log('Signal received: '+error.signal);
//    }
//    console.log('Child Process STDOUT: '+stdout);
//    console.log('Child Process STDERR: '+stderr);
//  });

function start(route){
  	function onRequest(request, response){
		var pathname = url.parse(request.url).pathname;
		console.log("Request for  " +pathname + "   received");
		// console.ls();
		
		route(pathname);
		console.log(pathname);
		if(pathname==="/shit")
		{
			var shit = require("./child");
		}
		else if (pathname === "/jeans")
		{
			var thing = require("./jean");
		}
		response.writeHead(200, {"Content-Type": "text/plain"});
		response.write("HELLO");
		response.end();
	}		
	http.createServer(onRequest).listen(8888);
	console.log("server start");				
}

//  ls.on('exit', function (code) {
//    console.log('Child process exited with exit code '+code);
//  });


exports.start = start;

//server.listen(8080);

console.log("server has started");
