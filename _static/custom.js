document.addEventListener("DOMContentLoaded", function() {
    // add links to nav boxes
    boxes = document.querySelectorAll(".toms-nav-container .box, .toms-nav-box");
    boxes.forEach(element => {
        element.addEventListener("click", function() {
            window.location.href = this.getAttribute("data-href");
        })
    });
});