import clsx from "clsx";
import Form from "@/components/Form";
import TextBox from "@/components/TextBox";
import classes from "./Feedback.module.css";

const Feedback = () => (
  <div className={clsx("col", classes.feedback)}>
    <Form onSubmit={console.log}>
      <TextBox label="Name" name="name" />
      <TextBox label="Username" name="username" />
      <TextBox label="Email" name="email" />

      <button type="submit">Submit</button>
    </Form>
  </div>
);

export default Feedback;
