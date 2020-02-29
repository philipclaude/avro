listener["initialize"] = function()
{
  listener.availableCommands = ["startserver","help","ls","clear"];
  listener.helpCommands = ["start the websocket server", "what this did","unsupported...","clear the terminal"];
}

isfile = function(text)
{
  if (text==".." || text=="../") return false;
  text = text.split(".");
  if (text.length>1) return true;
  return false;
}

listener["parseCommand"] = function(msg)
{
  var out = {text: [], color: ""};

  data = msg;
  data = data.split(" ");
  msg = data[0];

  // create a new line like a terminal when enter is hit
  if (msg=="") return out;

  if (msg=="send")
  {
  }
  else if (msg=="startserver")
  {
  }
  else if (msg=="killserver")
  {
    socket_killserver();
  }
  else if (msg=="help")
  {
    out = print_help(out);
  }
  else if (msg=="ls")
  {
    wv.socketUt.send("ls|");

    // the server will send a response back to wv and we can get the values of the directory
    var terminal = document.getElementById("brframe");
    terminal.listener.printDirectory();
  }
  else if (msg=="cd")
  {
    if (isfile(data[1]))
      alert("cannot cd into this!");
    else
      wv.socketUt.send("cd|"+data[1]);
  }
  else if (msg=="plot")
  {
    // need to do a lot of work to parse the rest of the command
    parsePlot(data);
  }
  else
  {
  }

  return out;
}

parsePlot = function(command)
{
  postMessage(command,"yellow",false,false);
  wv.socketUt.send("plot|"+command[1]); // for now...
}

listener["printDirectory"] = function()
{
  var drop = document.getElementById("directory");

  var text = "";
  for (var i=0;i<drop.options.length;i++)
    text = text +" "+drop.options[i].text;
  postMessage( text, "white" , false , false );

}

socket_sendmessage = function(msg)
{
}

socket_startserver = function()
{

}

socket_killserver = function()
{
}

print_help = function(out)
{
  for (i=0;i<listener.availableCommands.length;i++)
  {
    out.text[i] = listener.availableCommands[i]+": "+listener.helpCommands[i];
  }
  out.color = "yellow";
  return out;
}
