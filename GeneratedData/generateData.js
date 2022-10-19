function setup() {
    createCanvas(1000, 1000);
    noStroke();
    fill(0);
  }
  
  function draw() {
    // background(220);
    if (mouseIsPressed == true){
      ellipse(mouseX, mouseY, 3);
      console.log("${mouseX} ${mouseY}");
    }
  }