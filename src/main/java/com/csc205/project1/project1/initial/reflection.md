My AI of choice was chatGPT. Initially I used Claude but after the second prompt I hit a pay-wall.

The prompts:

Would you generate a class in Java that represents a line in 3D space. Utilize the Point3D class we created earlier as necessary. 
Include methods for line length, as well as a shortest distance between lines method and any other methods that make sense for a 3D line. Include comments for each method in the style of Spring Framework 
getting started guides explaining each method in detail. Also include a logger utilizing the Java util logging framework. Use log lines to print information and errors in each method that make sense for the method context. 
For now, just utilize the INFO, WARNING and SEVERE log levels. Document the object-oriented design patterns used in the class and explain how they demonstrate foundational principles for data structures and algorithms.

All of the errors I received were syntactically based. For example the in Point3D method return ChatGPT coded "X", "Y", "Z" and IntelliJ preffered "getX", "getY", "getZ". I discovered this by mousing over each error and following 
the IDE's suggested Context Action.

When uploading I had an issue with it being visible on the user end via the link provided by me. I solved this issue by reforking and cloning the work I had already pushed and committing & pushing it again. A new link was provided. 
I tested this link on a seperate web browser.
