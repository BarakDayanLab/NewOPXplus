import os.path

import cv2
import time


class BDWebCam:

    def __init__(self, device_index=0):
        self.device_index = device_index  # the camera device index (they may be more than one connected
        self.webcam = None
        pass

    def connect(self):
        # Try either VideoCapture(0) or VideoCapture(1) based on your camera availability
        self.webcam = cv2.VideoCapture(self.device_index)
        pass

    def disconnect(self):
        if self.webcam is not None:
            self.webcam.release()
        self.webcam = None
        pass

    def start_live_stream(self):

        while True:
            # Capture the video frame by frame
            ret, frame = self.webcam.read()

            # Display the resulting frame
            cv2.imshow('Live Stream (press q to quit)', frame)

            # the 'q' button is set as the
            # quitting button you may use any
            # desired button of your choice
            if cv2.waitKey(1) & 0xFF == ord('q'):
                break
        pass

    def snap(self):
        check, frame = self.webcam.read()
        return check, frame

    def timelapse(self, folder, delay=1, number_of_frames=-1):
        """
        Takes consequent images from webcam.

        delay - delay between frames (in N seconds)
        number_of_frames - how many images to take (-1 for endless)
        folder - where to save the images

        """
        count = 0
        while count < number_of_frames:
            # Capture the video frame by frame
            ret, frame = self.webcam.read()

            # Save image to folder
            current_date_time = time.strftime("%Y%m%d_%H%M%S")
            file_name = f'{current_date_time}_{count+1}.jpg'
            full_path = os.path.join(folder, file_name)
            cv2.imwrite(filename=full_path, img=frame)

            print(f'Took snapshot - {full_path}')
            count += 1
            time.sleep(delay)

        pass

    @staticmethod
    def test():

        webcam = BDWebCam()
        webcam.connect()

        mode = 'TIMELAPSE'

        if mode == 'LIVE':
            webcam.start_live_stream()
        elif mode == 'TIMELAPSE':
            webcam.timelapse(folder=r'C:\temp\timelapse', delay=10, number_of_frames=50)
        elif mode == 'SNAP':
            check, frame = webcam.snap()
            cv2.imshow("Image", frame)
            cv2.waitKey(0)

        webcam.disconnect()

        pass


        #cv2.imwrite(filename=r'<Your Directory>\saved_img.jpg', img=frame)

        webcam.release()
        pass

if __name__ == "__main__":

    BDWebCam.test()
    pass