var server = require("./server");
var router = require("./router");
// var shit = require("./child");

server.start(router.route);
